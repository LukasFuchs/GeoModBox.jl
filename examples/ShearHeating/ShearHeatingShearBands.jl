using Plots, ExtendableSparse
using GeoModBox
using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
using GeoModBox.AdvectionEquation.TwoD, GeoModBox.HeatEquation.TwoD
using GeoModBox.Tracers.TwoD
using Base.Threads
using Printf, LinearAlgebra
using TimerOutputs

function ShearHeatingShearBands()
    to          =   TimerOutput()
    @timeit to "Ini" begin
    # Define Initial Condition ========================================== #
    #   1) block
    Ini         =   (p=:ShearBandSetting,
                     V=:ShearBandPS,
                     ε = 5e-14,
    ) 
    radius      =   3.0e3           # [ m ]
    # ------------------------------------------------------------------- #
    # Plot Settings ===================================================== #
    Pl  =   (
        qinc    =   10,
        mainc   =   2,
        qsc     =   100*(60*60*24*365.25)
    )
    # ------------------------------------------------------------------- #
    # Geometry ========================================================== #
    M       =   Geometry(
        xmin    =   -35.0e3,    # [ m ]
        xmax    =   35.0e3,     # [ m ]
        ymin    =   0.0e3,      # [ m ]
        ymax    =   40.0e3,     # [ m ]
    )
    # -------------------------------------------------------------------- #
    # Grid =============================================================== #
    NC      =   (
        x   =   200, 
        y   =   100,
    )
    NV      =   (
        x   =   NC.x + 1,
        y   =   NC.y + 1,
    )
    Δ       =   GridSpacing(
        x   =   (M.xmax - M.xmin)/NC.x,
        y   =   (M.ymax - M.ymin)/NC.y,
    )
    x       =   (
        c   =   LinRange(M.xmin+Δ.x/2,M.xmax-Δ.x/2,NC.x),
        ce  =   LinRange(M.xmin - Δ.x/2.0, M.xmax + Δ.x/2.0, NC.x+2),
        v   =   LinRange(M.xmin,M.xmax,NV.x),
    )
    y       =   (
        c   =   LinRange(M.ymin+Δ.y/2,M.ymax-Δ.y/2,NC.y),
        ce  =   LinRange(M.ymin - Δ.x/2.0, M.ymax + Δ.x/2.0, NC.y+2),
        v   =   LinRange(M.ymin,M.ymax,NV.y),
    )
    x1      =   (
        c2d     =   x.c .+ 0*y.c',
        v2d     =   x.v .+ 0*y.v', 
        vx2d    =   x.v .+ 0*y.ce',
        vy2d    =   x.ce .+ 0*y.v',
    )
    x   =   merge(x,x1)
    y1      =   (
        c2d     =   0*x.c .+ y.c',
        v2d     =   0*x.v .+ y.v',
        vx2d    =   0*x.v .+ y.ce',
        vy2d    =   0*x.ce .+ y.v',
    )
    y   =   merge(y,y1)
    # -------------------------------------------------------------------- #
    # Physics ============================================================ #
    P   =   Physics(
        g       =   0.0,                #   Schwerebeschleunigung [m/s^2]
        ρ₀      =   2700.0,             #   Hintergunddichte [kg/m^3]
        k       =   2.5,                #   Thermische Leitfaehigkeit [ W/m/K ]
        cp      =   1050.0,             #   Heat capacity [ J/kg/K ]
        # α       =   2.0e-5,             #   Thermischer Expnasionskoef. [ K^-1 ]
        # Q₀      =   1.84e-08,           #   Waermeproduktionsrate pro Volumen [W/m^3]
        η₀      =   1e21,               #   Viskositaet [ Pa*s ]
        # ΔT      =   2500.0,             #   Temperaturdifferenz
        # Falls Ra < 0 gesetzt ist, dann wird Ra aus den obigen Parametern
        # berechnet. Falls Ra gegeben ist, dann wird die Referenzviskositaet so
        # angepasst, dass die Skalierungsparameter die gegebene Rayleigh-Zahl
        # ergeben.
        # Ra      =   1.0e6,              #   Rayleigh number
        # Ttop    =   273.15,             #   Temperatur an der Oberfläche [ K ]
    )
    # ------------------------------------------------------------------- #
    # g       =   0               #   Gravitational acceleration

    η₀      =   1.0e21          #   Reference Viscosity
    η₁      =   1.0e19          #   Block Viscosity
    ηᵣ      =   log10(η₁/η₀)
    η       =   [η₀,η₁]         #   Viscosity for phases

    ρ₀      =   2700.0             #   Background density
    ρ₁      =   2700.0             #   Block density
    ρ       =   [ρ₀,ρ₁] 

    phase   =   [0,1]
    # ------------------------------------------------------------------- #
    # Animation and Plot Settings ======================================= #
    path        =   string("./examples/ShearHeating/2D/Results/")
    save_fig    =   0
    anim        =   Plots.Animation(path, String[] )
    filename    =   string("ShearHeatingBands_",Ini.p,"_ηr_",round(ηᵣ),
                            "_tracers_DC")
    # ------------------------------------------------------------------- #
    # Allocation ======================================================== #
    D   =   DataFields(
        T       =   zeros(Float64,(NC...)),
        T0      =   zeros(Float64,(NC...)),
        T_ex    =   zeros(Float64,(NC.x+2,NC.y+2)),
        T_exo   =   zeros(Float64,(NC.x+2,NC.y+2)),
        vx      =   zeros(Float64,NV.x,NC.y+2),
        vy      =   zeros(Float64,NC.x+2,NV.y),
        Pt      =   zeros(Float64,NC...),
        p       =   zeros(Float64,NC...),
        p_ex    =   zeros(Float64,(NC.x+2,NC.y+2)),
        ρ       =   zeros(Float64,NC...),
        ρ_ex    =   zeros(Float64,(NC.x+2,NC.y+2)),
        vxc     =   zeros(Float64,NC...),
        vyc     =   zeros(Float64,NC...),
        vxco    =   zeros(Float64,(NC...)),
        vyco    =   zeros(Float64,(NC...)),
        vc      =   zeros(Float64,NC...),
        wt      =   zeros(Float64,(NC.x,NC.y)),
        wte     =   zeros(Float64,(NC.x+2,NC.y+2)),
        wtv     =   zeros(Float64,(NV.x,NV.y)),
        ηc      =   zeros(Float64,NC...),
        η_ex    =   zeros(Float64,(NC.x+2,NC.y+2)),
        ηv      =   zeros(Float64,NV...),
    )
    # Needed for the defect correction solution ---
    divV        =   zeros(Float64,NC...)
    ε           =   (
        xx      =   zeros(Float64,NC...), 
        yy      =   zeros(Float64,NC...), 
        xy      =   zeros(Float64,NV...),
    )
    τ           =   (
        xx      =   zeros(Float64,NC...), 
        yy      =   zeros(Float64,NC...), 
        xy      =   zeros(Float64,NV...),
    )
    # ------------------------------------------------------------------- #
    # Boundary Conditions =============================================== #
    VBC     =   (
        type    =   (E=:const,W=:const,S=:freeslip,N=:const),
        val     =   (E=zeros(NV.y),W=zeros(NV.y),S=zeros(NV.x),N=zeros(NV.x),
                    vxE=zeros(NC.y),vxW=zeros(NC.y),vyS=zeros(NC.x),vyN=zeros(NC.x)),
    )
    # ------------------------------------------------------------------- #
    # Initial Condition ================================================= #
    IniVelocity!(Ini.V,D,VBC,NV,Δ,M,x,y;Ini.ε)
    # ------------------------------------------------------------------- #
    # Time ============================================================== #
    T   =   ( 
        tmax    =   [0.0],  
        Δfac    =   1.0,    # Courant time factor, i.e. dtfac*dt_courant
        Δ       =   [0.0],
        time    =   [0.0,0.0],
    )
    T.tmax[1]   =   20.589 * 1e6 * (60*60*24*365.25)   # [ s ]
    nt          =   20
    # ------------------------------------------------------------------- #
    end
    # Tracer Advection ================================================== #
    @timeit to "Tracer Ini" begin
    nmx,nmy     =   3,3
    noise       =   0
    nmark       =   nmx*nmy*NC.x*NC.y
    Aparam      =   :phase
    MPC         =   (
            c       =   zeros(Float64,(NC.x,NC.y)),
            v       =   zeros(Float64,(NV.x,NV.y)),
            th      =   zeros(Float64,(nthreads(),NC.x,NC.y)),
            thv     =   zeros(Float64,(nthreads(),NV.x,NV.y)),
    )
    MAVG        = (
            PC_th   =   [similar(D.wte) for _ = 1:nthreads()],  # per thread
            PV_th   =   [similar(D.ηv) for _ = 1:nthreads()],   # per thread
            wte_th  =   [similar(D.wte) for _ = 1:nthreads()],  # per thread
            wtv_th  =   [similar(D.wtv) for _ = 1:nthreads()],  # per thread
    )
    Ma      =   IniTracer2D(Aparam,nmx,nmy,Δ,M,NC,noise,Ini.p,phase;
                        ellA=radius)
    # RK4 weights ---
    rkw     =   1.0/6.0*[1.0 2.0 2.0 1.0]   # for averaging
    rkv     =   1.0/2.0*[1.0 1.0 2.0 2.0]   # for time stepping
    # Count marker per cell ---
    CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,NV)
    # Interpolate from markers to cell ---
    Markers2Cells(Ma,nmark,MAVG.PC_th,D.ρ_ex,MAVG.wte_th,D.wte,x,y,Δ,Aparam,ρ)
    D.ρ     .=  D.ρ_ex[2:end-1,2:end-1]  
    Markers2Cells(Ma,nmark,MAVG.PC_th,D.p_ex,MAVG.wte_th,D.wte,x,y,Δ,Aparam,phase)
    D.p     .=  D.p_ex[2:end-1,2:end-1]
    Markers2Cells(Ma,nmark,MAVG.PC_th,D.η_ex,MAVG.wte_th,D.wte,x,y,Δ,Aparam,η)
    D.ηc    .=  D.η_ex[2:end-1,2:end-1]
    Markers2Vertices(Ma,nmark,MAVG.PV_th,D.ηv,MAVG.wtv_th,D.wtv,x,y,Δ,Aparam,η)
    end
    # System of Equations =============================================== #
    # Iterations --- 
    niter       =   50
    ϵ           =   1e-8
    # Numbering, without ghost nodes! ---
    off    = [  NV.x*NC.y,                          # vx
                NV.x*NC.y + NC.x*NV.y,              # vy
                NV.x*NC.y + NC.x*NV.y + NC.x*NC.y]  # Pt

    Num    =    (
        Vx  =   reshape(1:NV.x*NC.y, NV.x, NC.y), 
        Vy  =   reshape(off[1]+1:off[1]+NC.x*NV.y, NC.x, NV.y), 
        Pt  =   reshape(off[2]+1:off[2]+NC.x*NC.y,NC...),
                )
    δx      =   zeros(maximum(Num.Pt))
    F       =   zeros(maximum(Num.Pt))
    # Residuals ---
    Fm     =    (
        x       =   zeros(Float64,NV.x, NC.y), 
        y       =   zeros(Float64,NC.x, NV.y)
    )
    FPt     =   zeros(Float64,NC...)      
    # ------------------------------------------------------------------- #
    # Time Loop ========================================================= #
    @timeit to "Time Loop" begin
    for it = 1 # :nt
        # Update Time ---
        T.time[1]   =   T.time[2] 
        @printf("Time step: #%04d, Time [Myr]: %04e\n ",it,
                    T.time[1]/(60*60*24*365.25)/1.0e6)
        # Momentum Equation ===
        # Initial Residual ---------------------------------------------- #
        D.vx[2:end-1,:]    .=  0.0
        D.vy[:,1:end-1]    .=  0.0
        D.Pt               .=  0.0
        @timeit to "Solution Iteration" begin
        for iter = 1:niter
            @timeit to "Residual" begin
            Residuals2D!(D,VBC,ε,τ,divV,Δ,D.ηc,D.ηv,P.g,Fm,FPt)
            F[Num.Vx]   =   Fm.x[:]
            F[Num.Vy]   =   Fm.y[:]
            F[Num.Pt]   =   FPt[:]
            @printf("||R|| = %1.4e\n", norm(F)/length(F))
            norm(F)/length(F) < ϵ ? break : nothing
            end
            # Assemble Coefficients ===================================== #
            @timeit to "Assembly" begin
            K       =   Assembly(NC, NV, Δ, D.ηc, D.ηv, VBC, Num)
            end
            # ----------------------------------------------------------- #
            # Solution of the linear system ============================= #
            @timeit to "Solution" begin
            δx      =   - K \ F
            end
            # ----------------------------------------------------------- #
            # Update Unknown Variables ================================== #
            D.vx[:,2:end-1]     .+=  δx[Num.Vx]
            D.vy[2:end-1,:]     .+=  δx[Num.Vy]
            D.Pt                .+=  δx[Num.Pt]
        end
        end
        # --------------------------------------------------------------- #
        # Get the velocity on the centroids ---
        for i = 1:NC.x
            for j = 1:NC.y
                D.vxc[i,j]  = (D.vx[i,j+1] + D.vx[i+1,j+1])/2
                D.vyc[i,j]  = (D.vy[i+1,j] + D.vy[i+1,j+1])/2
            end
        end
        @. D.vc        = sqrt(D.vxc^2 + D.vyc^2)
        # ---
        @show(minimum(D.vc))
        @show(maximum(D.vc))
        # ---
        if T.time[2] >= T.tmax[1]
            it = nt
        end
        # ---
        if mod(it,2) == 0 || it == nt || it == 1
            # r = heatmap(x.c./1e3,y.c./1e3,D.ρ',color=:inferno,
            #             xlabel="x[km]",ylabel="y[km]",colorbar=false,
            #             title="ρ",
            #             aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),                             ylims=(M.ymin/1e3, M.ymax/1e3),
            #             layout=(3,1),subplot=1)
            p = scatter(Ma.x[1:Pl.mainc:end]./1e3,Ma.y[1:Pl.mainc:end]./1e3,
                        ms=1,ma=0.5,mc=Ma.phase[1:Pl.mainc:end],markerstrokewidth=0.0,
                        xlabel="x[km]",ylabel="y[km]",colorbar=false,
                        title="tracers",label="",
                        aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3), 
                        ylims=(M.ymin/1e3, M.ymax/1e3))
                        # layout=(3,1),subplot=2)
            q = heatmap(x.c./1e3,y.c./1e3,D.p',color=:inferno,
                        xlabel="x[km]",ylabel="y[km]",colorbar=false,
                        title="phase",
                        aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),                             
                        ylims=(M.ymin/1e3, M.ymax/1e3))
                        # layout=(3,1),subplot=3)
            r = heatmap(x.c./1e3,y.c./1e3,D.vc',
                        xlabel="x[km]",ylabel="y[km]",colorbar=false,
                        title="V_c",color=cgrad(:batlow),
                        aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),
                        ylims=(M.ymin/1e3, M.ymax/1e3))
                        # layout=(1,3),subplot=3)
            s = heatmap(x.v./1e3,y.ce./1e3,D.vx',
                        xlabel="x[km]",ylabel="y[km]",colorbar=true,
                        title="V_x",color=cgrad(:batlow),
                        aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),
                        ylims=(M.ymin/1e3, M.ymax/1e3))
                        # layout=(1,3),subplot=1)
            u = heatmap(x.ce./1e3,y.v./1e3,D.vy',
                        xlabel="x[km]",ylabel="y[km]",colorbar=true,
                        title="V_y",color=cgrad(:batlow),
                        aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),
                        ylims=(M.ymin/1e3, M.ymax/1e3))
                        # layout=(1,3),subplot=2)
            quiver!(r,x.c2d[1:Pl.qinc:end,1:Pl.qinc:end]./1e3,
                        y.c2d[1:Pl.qinc:end,1:Pl.qinc:end]./1e3,
                        quiver=(D.vxc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc,
                                D.vyc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc),        
                        la=0.5,color="white")
                        # layout=(1,3),subplot=3)
            t = heatmap(x.c./1e3,y.c./1e3,log10.(D.ηc'),color=reverse(cgrad(:roma)),
                        xlabel="x[km]",ylabel="y[km]",title="η_c",
                        # clims=(15,27),
                        aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3), 
                        ylims=(M.ymin/1e3, M.ymax/1e3),colorbar=true)
                        # layout=(2,2),subplot=3)
            if save_fig == 1
                Plots.frame(anim)
            elseif save_fig == 0
                display(p)
                display(q)
                display(r)
                display(t)
                display(s)
                display(u)
            end
        end
        if T.time[2] >= T.tmax[1]
            break
        end
        # # Calculate Time Stepping ---
        # T.Δ[1]      =   T.Δfac * minimum((Δ.x,Δ.y)) / 
        #                     (sqrt(maximum(abs.(D.vx))^2 + maximum(abs.(D.vy))^2))
        # @printf("\n")
        # # Calculate Time ---
        # T.time[2]   =   T.time[1] + T.Δ[1]
        # if T.time[2] > T.tmax[1] 
        #     T.Δ[1]      =   T.tmax[1] - T.time[1]
        #     T.time[2]   =   T.time[1] + T.Δ[1]
        # end
        # # Advection ===
        # @timeit to "Tracer Advection" begin
        # # Advect tracers ---
        # @printf("Running on %d thread(s)\n", nthreads())  
        # AdvectTracer2D(Ma,nmark,D,x,y,T.Δ[1],Δ,NC,rkw,rkv)
        # CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,NV)
        # @timeit to "Tracer Interpolation" begin
        # # Interpolate phase from tracers to grid ---
        # Markers2Cells(Ma,nmark,MAVG.PC_th,D.ρ_ex,MAVG.wte_th,D.wte,x,y,Δ,Aparam,ρ)
        # D.ρ     .=   D.ρ_ex[2:end-1,2:end-1]  
        # Markers2Cells(Ma,nmark,MAVG.PC_th,D.p_ex,MAVG.wte_th,D.wte,x,y,Δ,Aparam,phase)
        # D.p     .=  D.p_ex[2:end-1,2:end-1]
        # Markers2Cells(Ma,nmark,MAVG.PC_th,D.η_ex,MAVG.wte_th,D.wte,x,y,Δ,Aparam,η)
        # D.ηc    .=   D.η_ex[2:end-1,2:end-1]
        # Markers2Vertices(Ma,nmark,MAVG.PV_th,D.ηv,MAVG.wtv_th,D.wtv,x,y,Δ,Aparam,η)
        # end
        # end
    end # End Time Loop
    end
    # Save Animation ---
    if save_fig == 1
        # Write the frames to a GIF file
        Plots.gif(anim, string( path, filename, ".gif" ), fps = 15)
        foreach(rm, filter(startswith(string(path,"00")), readdir(path,join=true)))
    end
    display(to)
end

ShearHeatingShearBands()
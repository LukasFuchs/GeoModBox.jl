using Plots
using ExtendableSparse
using GeoModBox
using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
using GeoModBox.AdvectionEquation.TwoD
using GeoModBox.Tracers.TwoD
using Base.Threads
using Printf, LinearAlgebra

function RTI()
    save_fig    =   1
    # Define Initial Condition ========================================== #
    Ini         =   (p=:RTI,) 
    λ           =   3.0e3           #   Perturbation wavelength[ m ]
    δA          =   1500/15         #   Amplitude [ m ]
    # ------------------------------------------------------------------- #
    # Plot Settings ===================================================== #
    Pl  =   (
        qinc    =   4,
        mainc   =   2,
        qsc     =   100*(60*60*24*365.25)*3
    )
    # ------------------------------------------------------------------- #
    # Geometry ========================================================== #
    M       =   Geometry(
        ymin    =   -3.0e3,     # [ m ]
        ymax    =   0.0,
        xmin    =   0.0,
    )
    ar          =   Int64(2 * λ / (M.ymax-M.ymin))  #   aspect ratio
    M.xmax      =   (M.ymax-M.ymin)*ar
    # -------------------------------------------------------------------- #
    # Grid =============================================================== # 
    NC  =   (
        x   =   50*ar,
        y   =   50,
    )
    NV  =   (
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
    g       =   10.0                #   Gravitational acceleration [ m/s^2 ]
    # 0 - upper layer; 1 - lower layer
    η₀      =   1e19                #   Viscosity composition 0 [ Pa s ]
    η₁      =   1e13                #   Viscosity composition 1 [ Pa s]
    ηᵣ      =   log10(η₁/η₀)
    η       =   [η₀,η₁]             #   Viscosity for phases 

    ρ₀      =   3000.0              #   Density composition 0 [ kg/m^3 ]
    ρ₁      =   2900.0              #   Density composition 1 [ kg/m^3 ]
    ρ       =   [ρ₀,ρ₁]             #   Density for phases

    phase   =   [0,1]
    # -------------------------------------------------------------------- #
    # Animation and Plot Settings ======================================= #
    path        =   string("./examples/StokesEquation/2D/Results/")
    anim        =   Plots.Animation(path, String[] )
    filename    =   string(Ini.p,"_ηr_",round(ηᵣ),
                            "_tracers_DC")
    # ------------------------------------------------------------------- #
    # Allocation ======================================================== #
    D       =   (
        ρ       =   zeros(Float64,(NC...)),
        p       =   zeros(Float64,(NC...)),
        cp      =   zeros(Float64,(NC...)),
        vx      =   zeros(Float64,(NV.x,NV.y+1)),
        vy      =   zeros(Float64,(NV.x+1,NV.y)),    
        Pt      =   zeros(Float64,(NC...)),
        vxc     =   zeros(Float64,(NC...)),
        vyc     =   zeros(Float64,(NC...)),
        vc      =   zeros(Float64,(NC...)),
        wt      =   zeros(Float64,(NC.x,NC.y)),
        wtv     =   zeros(Float64,(NV.x,NV.y)),
        ηc      =   zeros(Float64,NC...),
        ηv      =   zeros(Float64,NV...),
    )
    # ------------------------------------------------------------------- #
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
        type    =   (E=:freeslip,W=:freeslip,S=:noslip,N=:noslip),
        val     =   (E=zeros(NV.y),W=zeros(NV.y),S=zeros(NV.x),N=zeros(NV.x)),
    )
    # ------------------------------------------------------------------- #
    # Time ============================================================== #
    T   =   TimeParameter(
        tmax    =   4500.0,         #   [ Ma ]
        Δfacc   =   1.0,            #   Courant time factor
        itmax   =   50,             #   Maximum iterations; 30000
    )
    T.tmax      =   T.tmax*1e6*T.year    #   [ s ]
    T.Δ         =   T.Δfacc * minimum((Δ.x,Δ.y)) / 
                        (sqrt(maximum(abs.(D.vx))^2 + maximum(abs.(D.vy))^2))

    Time        =   zeros(T.itmax)
    # ------------------------------------------------------------------- #
    # Tracer Advection ================================================== #
    nmx,nmy     =   5,5
    noise       =   0
    nmark       =   nmx*nmy*NC.x*NC.y
    Aparam      =   :phase
    MPC         =   (
        c       =   zeros(Float64,(NC.x,NC.y)),
        v       =   zeros(Float64,(NV.x,NV.y)),
        th      =   zeros(Float64,(nthreads(),NC.x,NC.y)),
        thv     =   zeros(Float64,(nthreads(),NV.x,NV.y)),
    )
    MPC1        = (
        PG_th   =   [similar(D.ρ) for _ = 1:nthreads()],    # per thread
        PV_th   =   [similar(D.ηv) for _ = 1:nthreads()],   # per thread
        wt_th   =   [similar(D.wt) for _ = 1:nthreads()],   # per thread
        wtv_th  =   [similar(D.wtv) for _ = 1:nthreads()],  # per thread
    )
    MPC     =   merge(MPC,MPC1)
    Ma      =   IniTracer2D(Aparam,nmx,nmy,Δ,M,NC,noise,Ini.p,phase;λ,δA)
    # RK4 weights ---
    rkw     =   1.0/6.0*[1.0 2.0 2.0 1.0]   # for averaging
    rkv     =   1.0/2.0*[1.0 1.0 2.0 2.0]   # for time stepping
    # Count marker per cell ---
    CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,NV,1)
    # Interpolate from markers to cell ---
    Markers2Cells(Ma,nmark,MPC.PG_th,D.ρ,MPC.wt_th,D.wt,x,y,Δ,Aparam,ρ)
    Markers2Cells(Ma,nmark,MPC.PG_th,D.ηc,MPC.wt_th,D.wt,x,y,Δ,Aparam,η)
    # Markers2Cells(Ma,nmark,MPC.PG_th,D.p,MPC.wt_th,D.wt,x,y,Δ,Aparam,phase)
    Markers2Vertices(Ma,nmark,MPC.PV_th,D.ηv,MPC.wtv_th,D.wtv,x,y,Δ,Aparam,η)
    # @. D.ηc     =   0.25 * (D.ηv[1:end-1,1:end-1] + 
    #                         D.ηv[2:end-0,1:end-1] + 
    #                         D.ηv[1:end-1,2:end-0] + 
    #                         D.ηv[2:end-0,2:end-0])
    # ------------------------------------------------------------------- #
    # System of Equations =============================================== #
    # Iterations
    niter   =   10
    ϵ       =   1e-8
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
    for it = 1:T.itmax
        # Update Time ---
        if it > 1
            Time[it]   =   Time[it-1] + T.Δ 
        end
        @printf("Time step: #%04d, Time [Myr]: %04e\n ",it,
                    Time[it]/(60*60*24*365.25)/1.0e6)
        # Momentum Equation ===
        # Initial Residual ---------------------------------------------- #
        D.vx    .=  0.0
        D.vy    .=  0.0
        D.Pt    .=  1.0
        for iter=1:niter
            Residuals2D!(D,VBC,ε,τ,divV,Δ,D.ηc,D.ηv,g,Fm,FPt)
            F[Num.Vx]   =   Fm.x[:]
            F[Num.Vy]   =   Fm.y[:]
            F[Num.Pt]   =   FPt[:]
            @printf("||R|| = %1.4e\n", norm(F)/length(F))
            norm(F)/length(F) < ϵ ? break : nothing
            # Assemble Coefficients ========================================= #
            K       =   Assembly(NC, NV, Δ, D.ηc, D.ηv, VBC, Num)
            # --------------------------------------------------------------- #
            # Solution of the linear system ================================= #
            δx      =   - K \ F
            # --------------------------------------------------------------- #
            # Update Unknown Variables ====================================== #
            D.vx[:,2:end-1]     .+=  δx[Num.Vx]
            D.vy[2:end-1,:]     .+=  δx[Num.Vy]
            D.Pt                .+=  δx[Num.Pt]
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
        if Time[it] >= T.tmax
            it = T.itmax
        end
        # ---
        if mod(it,2) == 0 || it == T.itmax || it == 1
            p = heatmap(x.c./1e3,y.c./1e3,D.ρ',color=:inferno,
                        xlabel="x[km]",ylabel="y[km]",colorbar=true,
                        title="ρ",
                        aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),                             
                        ylims=(M.ymin/1e3, M.ymax/1e3),
                        layout=(2,2),subplot=1)
            scatter!(p,Ma.x[1:Pl.mainc:end]./1e3,Ma.y[1:Pl.mainc:end]./1e3,
                        ms=1,ma=0.5,mc=Ma.phase[1:Pl.mainc:end],markerstrokewidth=0.0,
                        xlabel="x[km]",ylabel="y[km]",colorbar=true,
                        title="tracers",label="",
                        aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3), 
                        ylims=(M.ymin/1e3, M.ymax/1e3),
                        layout=(2,2),subplot=2)
            heatmap!(p,x.c./1e3,y.c./1e3,D.vc',
                        xlabel="x[km]",ylabel="y[km]",colorbar=true,
                        title="V_c",color=cgrad(:batlow),
                        aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),
                        ylims=(M.ymin/1e3, M.ymax/1e3),
                        layout=(2,2),subplot=4)
            quiver!(p,x.c2d[1:Pl.qinc:end,1:Pl.qinc:end]./1e3,
                        y.c2d[1:Pl.qinc:end,1:Pl.qinc:end]./1e3,
                        quiver=(D.vxc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc,
                                D.vyc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc),        
                        la=0.5,color="white",layout=(2,2),subplot=4)
            heatmap!(p,x.c./1e3,y.c./1e3,log10.(D.ηc'),color=reverse(cgrad(:roma)),
                        xlabel="x[km]",ylabel="y[km]",title="η_c",
                        clims=(15,27),
                        aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3), 
                        ylims=(M.ymin/1e3, M.ymax/1e3),colorbar=true,
                        layout=(2,2),subplot=3)
            if save_fig == 1
                Plots.frame(anim)
            elseif save_fig == 0
                display(p)
            end
        end
        if Time[it] >= T.tmax
            break
        end
        # Calculate Time Stepping ---
        T.Δ        =   T.Δfacc * minimum((Δ.x,Δ.y)) / 
                (sqrt(maximum(abs.(D.vx))^2 + maximum(abs.(D.vy))^2))
        if Time[it] > T.tmax
            T.Δ         =   T.tmax - Time[it-1]
            Time[it]    =   Time[it-1] + T.Δ
            it          =   T.itmax
        end
        # Advection ===
        # Advect tracers ---
        @printf("Running on %d thread(s)\n", nthreads())  
        AdvectTracer2D(Ma,nmark,D,x,y,T.Δ,Δ,NC,rkw,rkv,1)
        CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,NV,it)
        # Interpolate phase from tracers to grid ---
        Markers2Cells(Ma,nmark,MPC.PG_th,D.ρ,MPC.wt_th,D.wt,x,y,Δ,Aparam,ρ)
        Markers2Cells(Ma,nmark,MPC.PG_th,D.ηc,MPC.wt_th,D.wt,x,y,Δ,Aparam,η)
        # Markers2Cells(Ma,nmark,MPC.PG_th,D.p,MPC.wt_th,D.wt,x,y,Δ,Aparam,phase)
        Markers2Vertices(Ma,nmark,MPC.PV_th,D.ηv,MPC.wtv_th,D.wtv,x,y,Δ,Aparam,η)
        # @. D.ηc     =   0.25 * (D.ηv[1:end-1,1:end-1] + 
        #                     D.ηv[2:end-0,1:end-1] + 
        #                     D.ηv[1:end-1,2:end-0] + 
        #                     D.ηv[2:end-0,2:end-0])
    end # End Time Loop
    # Save Animation ---
    if save_fig == 1
        # Write the frames to a GIF file
        Plots.gif(anim, string( path, filename, ".gif" ), fps = 15)
        foreach(rm, filter(startswith(string(path,"00")), readdir(path,join=true)))
    end
end

RTI()
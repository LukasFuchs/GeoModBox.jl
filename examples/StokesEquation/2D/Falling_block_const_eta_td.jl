using Plots
using ExtendableSparse
using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
using GeoModBox.AdvectionEquation.TwoD
using GeoModBox.Tracers.TwoD
using Base.Threads
using Printf

function main()
    # Define Numerical Scheme =========================================== #
    # Advection ---
    #   1) upwind, 2) slf, 3) semilag, 4) tracers
    FD          =   (Method     = (Adv=:semilag,),)
    # ------------------------------------------------------------------- #
    # Define Initial Condition ========================================== #
    # Density --- 
    #   1) block
    Ini         =   (p=:block,) 
    # ------------------------------------------------------------------- #
    # Plot Settings ===================================================== #
    Pl  =   (
        qinc    =   5,
        qsc     =   100*(60*60*24*365.25)*1e0
    )
    # ------------------------------------------------------------------- #
    # Geometry ========================================================== #
    M       =   (
        xmin    =   0.0,
        xmax    =   500.0e3,    # [ m ]
        ymin    =   -500.0e3,   # [ m ]
        ymax    =   0.0,
    )
    # ------------------------------------------------------------------- #
    # Grid ============================================================== #
    NC      =   (
        x   =   50, 
        y   =   50,
    )
    NV      =   (
        x   =   NC.x + 1,
        y   =   NC.y + 1,
    )
    Δ       =   (
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
    # ------------------------------------------------------------------- #
    # Physics =========================================================== #
    g       =   9.81

    η₀      =   1.0e21

    ρ₀      =   3200.0          #   Background density
    ρ₁      =   3300.0          #   Block density
    ρ       =   [ρ₀,ρ₁] 

    phase   =   [0,1]
    # ------------------------------------------------------------------- #
    # Animationsettings ================================================== #
    path        =   string("./examples/StokesEquation/2D/Results/")
    anim        =   Plots.Animation(path, String[] )
    filename    =   string("Falling_",Ini.p,"_iso_",FD.Method.Adv)
    save_fig    =   1
    # -------------------------------------------------------------------- #
    # Allocation ======================================================== #
    D   =   (
        vx      =   zeros(Float64,NV.x,NC.y+2),
        vy      =   zeros(Float64,NC.x+2,NV.y),
        Pt      =   zeros(Float64,NC...),
        p       =   zeros(Int64,NC...),
        p_ex    =   zeros(Int64,NC.x+2,NC.y+2),
        ρ       =   zeros(Float64,NC...),
        ρ_ex    =   zeros(Float64,NC.x+2,NC.y+2),
        ρ_exo   =   zeros(Float64,NC.x+2,NC.y+2),
        vxc     =   zeros(Float64,NC...),
        vyc     =   zeros(Float64,NC...),
        vc      =   zeros(Float64,NC...),
        wt      =   zeros(Float64,(NC.x,NC.y)),
    )
    # ------------------------------------------------------------------- #
    # Boundary Conditions =============================================== #
    VBC     =   (
        type    =   (E=:freeslip,W=:freeslip,S=:freeslip,N=:freeslip),
        val     =   (E=zeros(NV.y),W=zeros(NV.y),S=zeros(NV.x),N=zeros(NV.x)),
    )
    # ------------------------------------------------------------------- #
    # Initial Condition ================================================= #
    # Phase ---
    IniPhase!(Ini.p,D,M,x,y,NC;phase)
    for i in eachindex(phase)
        D.ρ[D.p.==phase[i]] .= ρ[i]
    end
    D.ρ_ex[2:end-1,2:end-1]     .=  D.ρ
    D.ρ_ex[1,:]     .=   D.ρ_ex[2,:]
    D.ρ_ex[end,:]   .=   D.ρ_ex[end-1,:]
    D.ρ_ex[:,1]     .=   D.ρ_ex[:,2]
    D.ρ_ex[:,end]   .=   D.ρ_ex[:,end-1]
    # ------------------------------------------------------------------- #
    # Time ============================================================== #
    T   =   ( 
        tmax    =   [0.0],  
        Δfac    =   1.0,    # Courant time factor, i.e. dtfac*dt_courant
        Δ       =   [0.0],
        time    =   [0.0,0.0],
    )
    T.tmax[1]   =   9.886 * 1e6 * (60*60*24*365.25)   # [ s ]
    nt          =   9999
    # ------------------------------------------------------------------- #
    # Tracer Advection ================================================== #
    if FD.Method.Adv==:tracers 
        # Tracer Initialization ---
        nmx,nmy     =   3,3
        noise       =   1
        nmark       =   nmx*nmy*NC.x*NC.y
        Aparam      =   :phase
        MPC         =   (
            c               =   zeros(Float64,(NC.x,NC.y)),
            th              =   zeros(Float64,(nthreads(),NC.x,NC.y)),                
            min             =   zeros(Float64,nt),
            max             =   zeros(Float64,nt),
            mean            =   zeros(Float64,nt),
        )
        MPC1        = (
            PG_th   =   [similar(D.p) for _ = 1:nthreads()], # per thread
            wt_th   =   [similar(D.wt) for _ = 1:nthreads()], # per thread
        )
        MPC     =   merge(MPC,MPC1)
        Ma      =   IniTracer2D(Aparam,nmx,nmy,Δ,M,NC,noise)
        # RK4 weights ---
        rkw     =   1.0/6.0*[1.0 2.0 2.0 1.0]   # for averaging
        rkv     =   1.0/2.0*[1.0 1.0 2.0 2.0]   # for time stepping
        # Interpolate on centroids ---
        @threads for k = 1:nmark
            Ma.phase[k] =   FromCtoM(D.p_ex, k, Ma, x, y, Δ, NC)
        end
        # Count marker per cell ---
        CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,1)
    end
    # ------------------------------------------------------------------- #
    # System of Equations =============================================== #
    # Numbering, without ghost nodes! ---
    off    = [  NV.x*NC.y,                          # vx
                NV.x*NC.y + NC.x*NV.y,              # vy
                NV.x*NC.y + NC.x*NV.y + NC.x*NC.y]  # Pt

    Num    =    (
        Vx  =   reshape(1:NV.x*NC.y, NV.x, NC.y), 
        Vy  =   reshape(off[1]+1:off[1]+NC.x*NV.y, NC.x, NV.y), 
        Pt  =   reshape(off[2]+1:off[2]+NC.x*NC.y,NC...),
    )
    χ       =   zeros(maximum(Num.Pt))  #   Unknown Vector
    rhs     =   zeros(maximum(Num.Pt))  #   Right-hand Side
    # ------------------------------------------------------------------- #
    # Assemble Coefficients ============================================= #
    K       =   Assemblyc(NC, NV, Δ, η₀, VBC, Num)
    # ------------------------------------------------------------------- #
    # Time Loop ========================================================= #
    for it = 1:nt
        # Update Time ---
        T.time[1]   =   T.time[2] 
        @printf("Time step: #%04d, Time [Myr]: %04e\n ",it,T.time[1]/(60*60*24*365.25)/1.0e6)
        # Momentum Equation ===
        # Update RHS ---
        rhs     =   updaterhsc( NC, NV, Δ, η₀, D.ρ, g, VBC, Num, rhs )
        # Solve System of Equations ---
        χ       =   K \ rhs
        # Update Unknown Variables ---
        D.vx[:,2:end-1]     .=  χ[Num.Vx]
        D.vy[2:end-1,:]     .=  χ[Num.Vy]
        D.Pt                .=  χ[Num.Pt]
        # ===
        # Get the velocity on the centroids ---
        for i = 1:NC.x
            for j = 1:NC.y
                D.vxc[i,j]  = (D.vx[i,j+1] + D.vx[i+1,j+1])/2
                D.vyc[i,j]  = (D.vy[i+1,j] + D.vy[i+1,j+1])/2
            end
        end
        @. D.vc        = sqrt(D.vxc^2 + D.vyc^2)

        @show(maximum(D.vc))
        @show(minimum(D.Pt))
        @show(maximum(D.Pt))

        if T.time[2] >= T.tmax[1]
            it = nt
        end

        if mod(it,2) == 0 || it == nt || it == 1
            p = heatmap(x.c./1e3,y.c./1e3,D.ρ',color=:inferno,
                    xlabel="x[km]",ylabel="y[km]",colorbar=false,
                    title="Density",
                    aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3), 
                    ylims=(M.ymin/1e3, M.ymax/1e3),
                    layout=(2,2),subplot=1)
            quiver!(p,x.c2d[1:Pl.qinc:end,1:Pl.qinc:end]./1e3,
                        y.c2d[1:Pl.qinc:end,1:Pl.qinc:end]./1e3,
                        quiver=(D.vx[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc,
                                D.vyc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc),        
                        color="white",layout=(2,2),subplot=1)
            heatmap!(p,x.c./1e3,y.c./1e3,D.vxc',
                        xlabel="x[km]",ylabel="y[km]",colorbar=false,
                        title="V_x",
                        aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),
                        ylims=(M.ymin/1e3, M.ymax/1e3),
                        layout=(2,2),subplot=3)
            heatmap!(p,x.c./1e3,y.c./1e3,D.vyc',
                        xlabel="x[km]",ylabel="y[km]",colorbar=false,
                        title="V_y",
                        aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),
                        ylims=(M.ymin/1e3, M.ymax/1e3),
                        layout=(2,2),subplot=4)
            heatmap!(p,x.c./1e3,y.c./1e3,D.Pt',
                        xlabel="x[km]",ylabel="y[km]",colorbar=false,
                        title="P_t",
                        aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),
                        ylims=(M.ymin/1e3, M.ymax/1e3),
                        layout=(2,2),subplot=2)
            if save_fig == 1
                Plots.frame(anim)
            elseif save_fig == 0
                display(p)
            end
        end
        if T.time[2] >= T.tmax[1]
            break
        end
         # Calculate Time Stepping ---
        T.Δ[1]      =   T.Δfac * minimum((Δ.x,Δ.y)) / 
                            (sqrt(maximum(abs.(D.vx))^2 + maximum(abs.(D.vy))^2))
        
        @show T.Δ[1]
        @printf("\n")
        # Calculate Time ---
        T.time[2]   =   T.time[1] + T.Δ[1]
        if T.time[2] > T.tmax[1] 
            T.Δ[1]      =   T.tmax[1] - T.time[1]
            T.time[2]   =   T.time[1] + T.Δ[1]
        end
        # Advection ===
        if FD.Method.Adv==:upwind
            upwindc2D!(D.ρ,D.ρ_ex,D.vxc,D.vyc,NC,T.Δ[1],Δ.x,Δ.y)
        elseif FD.Method.Adv==:slf
            slfc2D!(D.ρ,D.ρ_ex,D.ρ_exo,D.vxc,D.vyc,NC,T.Δ[1],Δ.x,Δ.y)
        elseif FD.Method.Adv==:semilag
            semilagc2D!(D.ρ,D.ρ_ex,D.vxc,D.vyc,[],[],x,y,T.Δ[1])
        elseif FD.Method.Adv==:tracers
            # Advect tracers ---
            AdvectTracer2D(Ma,nmark,D,x,y,T.Δ[1],Δ,NC,rkw,rkv,1)
            CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,it)
          
            # Interpolate temperature from tracers to grid ---
            Markers2Cells(Ma,nmark,MPC.PG_th,D.p,MPC.wt_th,D.wt,x,y,Δ,Aparam)           
            D.p_ex[2:end-1,2:end-1]     .= D.p
        end
    end # End Time Loop
    # Save Animation ==================================================== #
    if save_fig == 1
        # Write the frames to a GIF file
        Plots.gif(anim, string( path, filename, ".gif" ), fps = 15)
        foreach(rm, filter(startswith(string(path,"00")), readdir(path,join=true)))
    end
    # ------------------------------------------------------------------- #
end

main()
using Plots, ExtendableSparse, LinearAlgebra
using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
using GeoModBox.AdvectionEquation.TwoD
using GeoModBox.Tracers.TwoD
using Base.Threads
using Printf, TimerOutputs, LaTeXStrings, Measures

# ======================================================================= #
# Helper function ======================================================= #
# ======================================================================= #
function UpdateRheology!(D,ρ,η,avg)
    # Density --- artihmetic averaging ---  
    @. D.ρ     =   ρ[1]*(1.0 - D.p) + ρ[2]*D.p
    D.ρ_ex[2:end-1,2:end-1]     .=  D.ρ
    D.ρ_ex[1,:]     .=  D.ρ_ex[2,:]
    D.ρ_ex[end,:]   .=  D.ρ_ex[end-1,:]
    D.ρ_ex[:,1]     .=  D.ρ_ex[:,2]
    D.ρ_ex[:,end]   .=  D.ρ_ex[:,end-1]
    # Viscosity - Centroids and vertices ---
    if avg ==     :arithmetic
        @. D.ηc =   (1.0 - D.p) * η[1] + D.p * η[2]
    elseif avg == :harmonic
        @. D.ηc =   1.0 / ( (1.0 - D.p) / η[1] + D.p / η[2] )
    elseif avg == :geometric
        @. D.ηc =   η[1]^(1.0 - D.p) * η[2]^D.p
    else
        error("Unknown viscosity averaging: $(avg)")
    end
    # --- Extended Centroids-
    D.η_ex[2:end-1,2:end-1]     .=  D.ηc
    D.η_ex[1,:]     .=  D.η_ex[2,:]
    D.η_ex[end,:]   .=  D.η_ex[end-1,:]
    D.η_ex[:,1]     .=  D.η_ex[:,2]
    D.η_ex[:,end]   .=  D.η_ex[:,end-1]
    # --- Vertices -
    if avg == :arithmetic
        @. D.ηv =
            0.25*(
                D.η_ex[1:end-1,1:end-1] + 
                D.η_ex[2:end  ,1:end-1] + 
                D.η_ex[1:end-1,2:end  ] + 
                D.η_ex[2:end  ,2:end  ]
            )
    elseif avg == :harmonic
        @. D.ηv =
            4.0/(
                1/D.η_ex[1:end-1,1:end-1] + 
                1/D.η_ex[2:end  ,1:end-1] + 
                1/D.η_ex[1:end-1,2:end  ] + 
                1/D.η_ex[2:end  ,2:end  ]
            )
    elseif avg == :geometric
        @. D.ηv =
            exp(0.25*(
                log(D.η_ex[1:end-1,1:end-1]) + 
                log(D.η_ex[2:end  ,1:end-1]) + 
                log(D.η_ex[1:end-1,2:end  ]) + 
                log(D.η_ex[2:end  ,2:end  ])
            ))
    else
        error("Unknown viscosity averaging: $(avg)")
    end
    
    return D

end
# ======================================================================= #

# ======================================================================= #
# Start Main Function =================================================== #
# ======================================================================= #
function FallingBlockBenchmark(td)
    to      =   TimerOutput()
    @timeit to "Ini" begin
    # Benchmark parameter =============================================== #
    ηᵣ      =   LinRange(-6.0,6.0,13)       #   Viscosity ratio
    sv      =   zeros(length(ηᵣ))           #   Sinking Velocity
    tmax    =   [7.115094, 7.114844, 7.256534, 7.377311, 7.738412, 
                    7.673613, 9.886, 15.446, 19.623, 20.569, 20.569,
                    20.569, 20.589]
    # ------------------------------------------------------------------- #
    # Define Numerical Scheme =========================================== #
    # Advection ---
    #   1) upwind, 2) slf, 3) semilag, 4) tracers
    # Advection methods:
    #
    # :tracers
    #     Recommended for this benchmark. Preserves the sharp material
    #     interface and remains robust for large viscosity contrasts.
    #
    # :semilag
    #     Produces a comparatively coherent interface but requires bounded
    #     phase values and careful departure-point boundary treatment.
    #
    # :upwind
    #     Robust and bounded, but strongly diffuses the material interface.
    #
    # :slf
    #     Not recommended for discontinuous phase fields because dispersive
    #     oscillations and the leapfrog computational mode destabilize the
    #     interface.
    FD          =   (Method     = (Adv=:tracers,),)
    # Arithmetic averaging is used for this benchmark because harmonic and
    # geometric mixing produce strongly reduced mixed-cell viscosities for
    # weak inclusions. In combination with semi-Lagrangian phase transport,
    # this leads to large local velocities, very small CFL timesteps, and a
    # substantial increase in the required number of Stokes solves.
    #
    # This is an implementation-specific choice and not a general statement
    # that arithmetic viscosity averaging is optimal for all multiphase
    # Stokes problems.
    #
    # Eulerian methods:
    #   `avg` controls viscosity reconstruction from the transported phase
    #   ratio and interpolation from cell centroids to vertices.
    #
    # Tracer method:
    #   `avgm` controls direct averaging of marker viscosities onto the grid.
    avg         =   :arithmetic
    avgm        =   :arith
    # ------------------------------------------------------------------- #
    # Define Initial Condition ========================================== #
    # Density --- 
    #   1) block
    Ini         =   (p=:block,) 
    # ------------------------------------------------------------------- #
    # Animation and Plot Settings ======================================= #
    path        =   string("./examples/StokesEquation/2D/Results/")
    save_fig    =   1
    p2          =   plot(0,0,layout=(2,3))
    count       =   Int64(0)
    panel       =   ["(a)","(b)","(c)","(d)","(e)","(f)",]
    # ------------------------------------------------------------------- #
    # Plot Settings ===================================================== #
    Pl  =   (
        qinc    =   1,
        qsc     =   100*(60*60*24*365.25)*5e1
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
        ce  =   LinRange(M.ymin - Δ.y/2.0, M.ymax + Δ.y/2.0, NC.y+2),
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
    g       =   9.81                #   Gravitational Acceleration

    η₀      =   1.0e21              #   Background Viscosity
    
    ρ₀      =   3200.0              #   Background density
    ρ₁      =   3300.0              #   Block density
    ρ       =   [ρ₀,ρ₁]             #   Density for phases

    phase   =   [0,1]               #   Phase ID
    # ------------------------------------------------------------------- #
    # Boundary Conditions =============================================== #
    VBC     =   (
        type    =   (E=:freeslip,W=:freeslip,S=:freeslip,N=:freeslip),
        val     =   (E=zeros(NV.y),W=zeros(NV.y),S=zeros(NV.x),N=zeros(NV.x),
                    vxE=zeros(NC.y),vxW=zeros(NC.y),vyS=zeros(NC.x),vyN=zeros(NC.x)),
    )
    # ------------------------------------------------------------------- #
    end
    @timeit to "η loop" begin
    for mn in eachindex(ηᵣ)     #   Loop over ηᵣ
        anim        =   Plots.Animation(path, String[] )
        filename    =   string("Falling_",Ini.p,"_ηr_",round(ηᵣ[mn]),
                            "_",FD.Method.Adv,"_dc_",avg)
        # --------------------------------------------------------------- #
        # Physics ======================================================= #
        η₁      =   η₀ * 10^(ηᵣ[mn])    #   Block Viscosity
        η       =   [η₀,η₁]             #   Viscosity for phases
        @show η
        # --------------------------------------------------------------- #
        # Allocation ==================================================== #
        D   =   (
            vx      =   zeros(Float64,NV.x,NC.y+2),
            vy      =   zeros(Float64,NC.x+2,NV.y),
            Pt      =   zeros(Float64,NC...),
            p       =   zeros(Float64,NC...),
            p_ex    =   zeros(Float64,NC.x+2,NC.y+2),
            p_exo   =   zeros(Float64,NC.x+2,NC.y+2),
            ρ       =   zeros(Float64,NC...),
            ρ_ex    =   zeros(Float64,NC.x+2,NC.y+2),
            ρ_exo   =   zeros(Float64,NC.x+2,NC.y+2),
            vxc     =   zeros(Float64,NC...),
            vyc     =   zeros(Float64,NC...),
            vxco    =   zeros(Float64,NC...),
            vyco    =   zeros(Float64,NC...),
            vc      =   zeros(Float64,NC...),
            wt      =   zeros(Float64,(NC.x,NC.y)),
            wte     =   zeros(Float64,(NC.x+2,NC.y+2)),
            wtv     =   zeros(Float64,(NV.x,NV.y)),
            ηc      =   zeros(Float64,NC...),
            ηv      =   zeros(Float64,NV...),
            η_ex    =   zeros(Float64,NC.x+2,NC.y+2),
            η_exo   =   zeros(Float64,NC.x+2,NC.y+2),
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
        # --------------------------------------------------------------- #
        # Time ========================================================== #
        if FD.Method.Adv==:semilag
            Δfac    =   0.25
        else
            Δfac    =   0.9
        end
        T   =   ( 
            tmax    =   [0.0],  
            Δfac    =   Δfac,    # Courant time factor, i.e. dtfac*dt_courant
            Δ       =   [0.0],
            time    =   [0.0,0.0],
        )
        T.tmax[1]   =   tmax[mn] * 1e6 * (60*60*24*365.25)   # [ s ] 
        if td == 0
            nt  =   1
        else
            nt  =   9999
        end
        # --------------------------------------------------------------- #
        # Tracer Advection ============================================== #
        if FD.Method.Adv==:tracers 
            @timeit to "Tracer Ini" begin
            # Tracer Initialization ---
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
            Ma      =   IniTracer2D(Aparam,nmx,nmy,Δ,M,NC,noise,Ini.p,phase)
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
            Markers2Cells(Ma,nmark,MAVG.PC_th,D.η_ex,MAVG.wte_th,D.wte,x,y,Δ,Aparam,η;avgm=avgm)
            D.ηc    .=  D.η_ex[2:end-1,2:end-1]
            Markers2Vertices(Ma,nmark,MAVG.PV_th,D.ηv,MAVG.wtv_th,D.wtv,x,y,Δ,Aparam,η;avgm=avgm)
            end
        else
            @timeit to "Phase Ini" begin
            # ----------------------------------------------------------- #
            # Initial Condition ========================================= #
            IniPhase!(Ini.p,D,M,x,y,NC;phase)
            UpdateRheology!(D,ρ,η,avg)
            end
        end
        # --------------------------------------------------------------- #
        # System of Equations =========================================== #
        # Iterations --- 
        niter      =   50
        atol       =   1e-8        #   Absolute tolerance
        rtol       =   1e-5        #   Relative tolerance; r = atolM0/atolM
        RM         =   0.0         #   Initialize absolute residual    
        RMrel      =   0.0         #   Initialize relative residual 
        # Numbering, without ghost nodes! ---
        off    = [  NV.x*NC.y,                          # vx
                    NV.x*NC.y + NC.x*NV.y,              # vy
                    NV.x*NC.y + NC.x*NV.y + NC.x*NC.y]  # Pt

        Num    =    (
            Vx  =   reshape(1:NV.x*NC.y, NV.x, NC.y), 
            Vy  =   reshape(off[1]+1:off[1]+NC.x*NV.y, NC.x, NV.y), 
            Pt  =   reshape(off[2]+1:off[2]+NC.x*NC.y,NC...),
        )
        ndof    =   maximum(Num.Pt)        
        K       =   ExtendableSparseMatrix(ndof,ndof)      
        δx      =   zeros(maximum(Num.Pt))
        F       =   zeros(maximum(Num.Pt))
            # Residuals ---
        Fm     =    (
            x       =   zeros(Float64,NV.x, NC.y), 
            y       =   zeros(Float64,NC.x, NV.y)
        )
        FPt     =   zeros(Float64,NC...)      
        # --------------------------------------------------------------- #
        @timeit to "Time Loop" begin
        # Time Loop ===================================================== #
        for it = 1:nt
            R0      =   0.0
            # Update Time ---
            T.time[1]   =   T.time[2] 
            @printf("Time step: #%04d, Time [Myr]: %04e\n ",it,
                        T.time[1]/(60*60*24*365.25)/1.0e6)
            # Initial Residual ---------------------------------------------- #
            D.vx    .=  0.0
            D.vy    .=  0.0
            D.Pt    .=  0.0
            K       =   Assembly(NC, NV, Δ, D.ηc, D.ηv, VBC, Num)
            Kfac    =   lu(K.cscmatrix)
            for iter = 1:niter
                Residuals2D!(D,VBC,ε,τ,divV,Δ,D.ηc,D.ηv,g,Fm,FPt)
                F[Num.Vx]   .=   Fm.x
                F[Num.Vy]   .=   Fm.y
                F[Num.Pt]   .=   FPt
                RM          =   norm(F)/length(F)
                if iter == 1
                    R0 = RM
                end
                RMrel       =   RM/R0
                # if verbose_step
                    @printf("   MCE %2d: ||R|| = %1.4e, ||R||/||R₀|| = %1.4e\n",iter,RM,RMrel)
                # end
                (RM < atol || RM/R0 < rtol) && break
                δx      .=   - (Kfac \ F)
                # ----------------------------------------------------------- #
                # Update Unknown Variables ================================== #
                D.vx[:,2:end-1]     .+=  δx[Num.Vx]
                D.vy[2:end-1,:]     .+=  δx[Num.Vy]
                D.Pt                .+=  δx[Num.Pt]
            end
            # ======
            # Get the velocity on the centroids ---
            for i = 1:NC.x
                for j = 1:NC.y
                    D.vxc[i,j]  = (D.vx[i,j+1] + D.vx[i+1,j+1])/2
                    D.vyc[i,j]  = (D.vy[i+1,j] + D.vy[i+1,j+1])/2
                end
            end
            if FD.Method.Adv==:semilag
                if it == 1
                    @. D.vxco   =   D.vxc
                    @. D.vyco   =   D.vyc
                end
            end
            @. D.vc        = sqrt(D.vxc^2 + D.vyc^2)
            # ---
            if it == 1
                sv[mn]  =   maximum(D.vc)
            end
            # ---
            if T.time[2] >= T.tmax[1]
                it = nt
            end
            # ---
            if mod(it,2) == 0 || it == nt || it == 1
                if FD.Method.Adv==:tracers
                    p = heatmap(x.c./1e3,y.c./1e3,D.p',color=:inferno,
                                xlabel="x[km]",ylabel="y[km]",colorbar=false,
                                title="Phase_c",
                                aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3), 
                                ylims=(M.ymin/1e3, M.ymax/1e3),
                                layout=(2,2),subplot=1)
                else
                    p = heatmap(x.v./1e3,y.v./1e3,log10.(abs.(D.ηv')),color=reverse(cgrad(:roma)),
                                xlabel="x[km]",ylabel="y[km]",title="η_v",
                                clims=(15,27),
                                aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3), 
                                ylims=(M.ymin/1e3, M.ymax/1e3),colorbar=true,
                                layout=(2,2),subplot=1)
                end
                if FD.Method.Adv==:tracers
                    scatter!(p,Ma.x[1:Pl.qinc:end]./1e3,Ma.y[1:Pl.qinc:end]./1e3,
                                ms=1,ma=0.5,mc=Ma.phase[1:Pl.qinc:end],markerstrokewidth=0.0,
                                xlabel="x[km]",ylabel="y[km]",colorbar=false,
                                title="tracers",label="",
                                aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3), 
                                ylims=(M.ymin/1e3, M.ymax/1e3),
                                layout=(2,2),subplot=2)
                else
                    heatmap!(p,x.c./1e3,y.c./1e3,D.ρ',color=:inferno,
                                xlabel="x[km]",ylabel="y[km]",colorbar=false,
                                title="Density",
                                aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3), 
                                ylims=(M.ymin/1e3, M.ymax/1e3),
                                layout=(2,2),subplot=2)
                end
                heatmap!(p,x.c./1e3,y.c./1e3,D.vc',
                            xlabel="x[km]",ylabel="y[km]",colorbar=false,
                            title="V_c",color=cgrad(:batlow),
                            aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),
                            ylims=(M.ymin/1e3, M.ymax/1e3),
                            layout=(2,2),subplot=4)
                quiver!(p,x.c2d[1:Pl.qinc:end,1:Pl.qinc:end]./1e3,
                            y.c2d[1:Pl.qinc:end,1:Pl.qinc:end]./1e3,
                            quiver=(D.vxc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc,
                                    D.vyc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc),        
                            la=0.5,color="white",layout=(1,3),subplot=3)
                heatmap!(p,x.c./1e3,y.c./1e3,log10.(abs.(D.ηc)'),color=reverse(cgrad(:roma)),
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
            if T.time[2] >= T.tmax[1]
                break
            end
            # Calculate Time Stepping ---
            T.Δ[1]      =   T.Δfac * minimum((Δ.x,Δ.y)) / 
                                (sqrt(maximum(abs.(D.vx))^2 + maximum(abs.(D.vy))^2))
            # ---
            @printf("\n")
            # Calculate Time ---
            T.time[2]   =   T.time[1] + T.Δ[1]
            if T.time[2] > T.tmax[1] 
                T.Δ[1]      =   T.tmax[1] - T.time[1]
                T.time[2]   =   T.time[1] + T.Δ[1]
            end
            # Advection ======
            if FD.Method.Adv==:upwind
                upwindc2D!(D.p,D.p_ex,D.vxc,D.vyc,NC,T.Δ[1],Δ.x,Δ.y)
            elseif FD.Method.Adv==:slf
                slfc2D!(D.p,D.p_ex,D.p_exo,D.vxc,D.vyc,NC,T.Δ[1],Δ.x,Δ.y)
            elseif FD.Method.Adv==:semilag
                semilagc2D!(D.p,D.p_ex,D.vxc,D.vyc,D.vxco,D.vyco,x,y,T.Δ[1])
            elseif FD.Method.Adv==:tracers
                # Advect tracers ---
                @printf("Running on %d thread(s)\n", nthreads())  
                AdvectTracer2D(Ma,nmark,D,x,y,T.Δ[1],Δ,NC,rkw,rkv)
                CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,NV)
                # Interpolate phase from tracers to grid ---
                Markers2Cells(Ma,nmark,MAVG.PC_th,D.ρ_ex,MAVG.wte_th,D.wte,x,y,Δ,Aparam,ρ)
                D.ρ     .=   D.ρ_ex[2:end-1,2:end-1]  
                Markers2Cells(Ma,nmark,MAVG.PC_th,D.p_ex,MAVG.wte_th,D.wte,x,y,Δ,Aparam,phase)
                D.p     .=  D.p_ex[2:end-1,2:end-1]
                Markers2Cells(Ma,nmark,MAVG.PC_th,D.η_ex,MAVG.wte_th,D.wte,x,y,Δ,Aparam,η;avgm=avgm)
                D.ηc    .=   D.η_ex[2:end-1,2:end-1]
                Markers2Vertices(Ma,nmark,MAVG.PV_th,D.ηv,MAVG.wtv_th,D.wtv,x,y,Δ,Aparam,η;avgm=avgm)
            end
            if FD.Method.Adv!=:tracers
                @. D.p      = clamp(D.p, 0.0, 1.0)
                @. D.p_ex   = clamp(D.p_ex, 0.0, 1.0)
                @. D.p_exo  = clamp(D.p_exo, 0.0, 1.0)
                UpdateRheology!(D,ρ,η,avg)
            end
            @printf("\n")
            if FD.Method.Adv==:semilag
                @. D.vxco   =   D.vxc
                @. D.vyco   =   D.vyc
            end
        end     # End Time Loop
        end
        if ηᵣ[mn] == 0.0 || ηᵣ[mn] == 1.0 || ηᵣ[mn] == 2.0 || 
                        ηᵣ[mn] == 3.0 || ηᵣ[mn] == 4.0 || ηᵣ[mn] == 6.0
            count = count + 1
            if FD.Method.Adv==:tracers
                if ηᵣ[mn] == 3.0 || ηᵣ[mn] ==  4.0 || ηᵣ[mn] == 6.0
                    xlab    =   L"x[km]"
                    xf      =   :auto
                else
                    xlab    =   ""
                    xf      =   _ -> ""
                end
                if ηᵣ[mn] == 0.0 || ηᵣ[mn] == 3.0
                    ylab    =   L"y[km]"
                    yf      =   :auto
                else
                    ylab    =   ""
                    yf      =   _ -> ""
                end
                p2  =   scatter!(p2,Ma.x[1:Pl.qinc:end]./1e3,Ma.y[1:Pl.qinc:end]./1e3,
                    ms=2,ma=0.5,mc=Ma.phase[1:Pl.qinc:end],markerstrokewidth=0.0,
                    xlabel=xlab,ylabel= ylab,colorbar=false,
                    title           = latexstring(panel[count],"\\quad\\log_{10}(\\eta_r)=",string(ηᵣ[mn])),
                    titlefontsize   = 20,
                    label="", xformatter = xf,yformatter = yf,
                    aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3), 
                    ylims=(M.ymin/1e3, M.ymax/1e3),
                    layout=(2,3),subplot=count,
                    framestyle      = :box,
                    guidefontsize   = 20,
                    tickfontsize    = 16,
                    size            = (1200, 800),
                    right_margin    = 8mm,
                    left_margin     = 8mm,
                    dpi             = 300)
            else
                p2 = heatmap!(p2,x.c./1e3,y.c./1e3,D.p',color=:inferno,
                    xlabel="x[km]",ylabel="y[km]",colorbar=false,
                    title="Phase_c",
                    aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3), 
                    ylims=(M.ymin/1e3, M.ymax/1e3),
                    layout=(2,3),subplot=count)
            end
        end
        # Save Animation ---
        if save_fig == 1
            if td == 1
                # Write the frames to a GIF file
                Plots.gif(anim, string( path, filename, ".gif" ), fps = 15)
                foreach(rm, filter(startswith(string(path,"00")), readdir(path,join=true)))
            end
        end
    end # End ηᵣ Loop
    end
    if td == 0
        q = scatter(ηᵣ,sv.*(100.0*365.25*24*60*60),
                        ylabel          = L"v_\mathrm{block} [cm/a]",
                        xlabel          = L"log_{10}\left(\eta_\mathrm{b}/\eta_\mathrm{m}\right)",
                        title           = "Sinking Velocity",
                        label           = false,
                        marker          = :circle,
                        markersize      = 7,
                        markercolor     = :black,
                        markerstrokecolor = :black,
                        markerstrokewidth = 0.8,
                        framestyle      = :box,
                        grid            = false,
                        guidefontsize   = 20,
                        tickfontsize    = 16,
                        legendfontsize  = 14,
                        xlims           = (-6, 6),
                        ylims           = (1.0, 5.0),
                        xticks          = -6:2:6,
                        size            = (900, 650),
                        left_margin     = 5mm,
                        right_margin    = 5mm,
                        bottom_margin   = 5mm,
                        top_margin      = 5mm,
                        dpi             = 300,)
        plot!(
            q,
            ηᵣ,
            sv .*(100.0*365.25*24*60*60),
            color      = :black,
            linewidth  = 2.0,
            linestyle  = :solid,
            label      = false,
        )
        if save_fig == 1
            savefig(q,string("./examples/StokesEquation/2D/Results/FallingBlock_SinkingVeloc",
                                "_",FD.Method.Adv,"_dc_",avg,".png"))
            foreach(rm, filter(startswith(string(path,"00")), readdir(path,join=true)))
        else
            display(q)
        end
    else
        if save_fig == -1 ||save_fig == 1
            savefig(p2,string("./examples/StokesEquation/2D/Results/FallingBlock_FinalStage",
                                "_",FD.Method.Adv,"_dc_",avg,".png"))
        else
            display(p2)
        end
    end
    display(to)
end
# ======================================================================= #
# Define if the problem is time-dependent (1) or if you want to have the  #
# steady state (0) solution.                                              #
td  =   1
# ---
FallingBlockBenchmark(td)
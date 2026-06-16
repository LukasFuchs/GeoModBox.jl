using Plots, ExtendableSparse
using GeoModBox
using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
using GeoModBox.AdvectionEquation.TwoD, GeoModBox.HeatEquation.TwoD
using GeoModBox.Tracers.TwoD
using Base.Threads
using Printf, LinearAlgebra
using TimerOutputs
function GetSecondInvariants!(ε,τ)
    # Get shear strain rate and stress on centroids ---
    @. ε.xyc    =   (ε.xy[1:end-1,1:end-1] + ε.xy[2:end-0,1:end-1] + 
                    ε.xy[1:end-1,2:end-0] + ε.xy[2:end-0,2:end-0])/4.0
    @. τ.xyc    =   (τ.xy[1:end-1,1:end-1] + τ.xy[2:end-0,1:end-1] + 
                    τ.xy[1:end-1,2:end-0] + τ.xy[2:end-0,2:end-0])/4.0
    # Get invariants ---
    @. ε.II     =   sqrt(0.5*(ε.xx^2 + ε.yy^2) + ε.xyc^2)
    @. τ.II     =   sqrt(0.5*(τ.xx^2 + τ.yy^2) + τ.xyc^2)
end

function ShearHeatingShearBands()
    to          =   TimerOutput()
    @timeit to "Ini" begin
    # Define Initial Condition ========================================== #
    #   1) block
    Ini         =   (T=:const,
                     p=:ShearBandSetting,
                     V=:ShearBandPS,
                     ε = 5e-14,
    ) 
    radius      =   3.0e3           # [ m ]
    strain      =   0
    # ------------------------------------------------------------------- #
    # Define numerical methods ========================================== #
    FD  =   (Method = (
            Diff = :dc,
            Adv  = :slf), 
    )
    # if Diff =: dc; θ-rule
    #       C = 0   -> implicit
    #       C = 0.5 -> CN discretization
    #       C = 1.0 -> explicit
    C   =   0.5
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
    # -------------------------------------------------------------------- #
    # Physics ============================================================ #
    P   =   Physics(
        g       =   0.0,            #   Gravitational acceleration [m/s^2]
        ρ₀      =   2700.0,         #   Reference density [kg/m^3]
        k       =   2.5,            #   Thermal conductivity [ W/m/K ]
        cp      =   1050.0,         #   Heat capacity [ J/kg/K ]
        η₀      =   1e21,           #   Reference viscosity [ Pa*s ]
    )
    # P.κ         =   0.0
    # ------------------------------------------------------------------- #
    # g       =   0               #   Gravitational acceleration

    η₀      =   1.0e20          #   Reference Viscosity
    η₁      =   1.0e19          #   Block Viscosity
    ηᵣ      =   log10(η₁/η₀)
    η       =   [η₀,η₁]         #   Viscosity for phases

    ρ₀      =   2700.0          #   Background density
    ρ₁      =   2700.0          #   Block density
    ρ       =   [ρ₀,ρ₁] 

    phase   =   [0,1]
    # Define rheology paramters ========================================= #
    # ------------------------------------------------------------------- #
    # ------------------------------------------------------------------- #
    # Animation and Plot Settings ======================================= #
    path        =   string("./examples/ShearHeating/2D/Results/")
    save_fig    =   1
    anim        =   Plots.Animation(path, String[] )
    filename    =   string("ShearHeatingBands")
    # ------------------------------------------------------------------- #
    # Allocation ======================================================== #
    D   =   DataFields(
        T       =   zeros(Float64,(NC...)),
        T0      =   zeros(Float64,(NC...)),
        T_ex    =   zeros(Float64,(NC.x+2,NC.y+2)),
        T_ex0   =   zeros(Float64,(NC.x+2,NC.y+2)),
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
    # Needed for Euler advection schemes ---
    vx_bg_c     =   zeros(Float64,NC...)
    vy_bg_c     =   zeros(Float64,NC...)
    vxc_res     =   zeros(Float64,NC...)
    vyc_res     =   zeros(Float64,NC...)
    vxco_res    =   zeros(Float64,NC...)
    vyco_res    =   zeros(Float64,NC...)
    # Needed for the defect correction solution ---
    divV        =   zeros(Float64,NC...)
    ε           =   (
        xx      =   zeros(Float64,NC...), 
        yy      =   zeros(Float64,NC...), 
        xy      =   zeros(Float64,NV...),
        xyc     =   zeros(Float64,NC...),
        II      =   zeros(Float64,NC...),
    )
    τ           =   (
        xx      =   zeros(Float64,NC...), 
        yy      =   zeros(Float64,NC...), 
        xy      =   zeros(Float64,NV...),
        xyc     =   zeros(Float64,NC...),
        II      =   zeros(Float64,NC...),
    )
    Hₛ      =   zeros(Float64,NC...)
    # ------------------------------------------------------------------- #
    # Boundary Conditions =============================================== #
    VBC     =   (
        type    =   (E=:const,W=:const,S=:freeslip,N=:const),
        val     =   (E=zeros(NV.y),W=zeros(NV.y),S=zeros(NV.x),N=zeros(NV.x),
                    vxE=zeros(NC.y),vxW=zeros(NC.y),vyS=zeros(NC.x),vyN=zeros(NC.x)),
    )
    TBC     =   (
        type    =   (E=:Neumann,W=:Neumann,S=:Neumann,N=:Neumann),
        val     =   (E=zeros(NC.y),W=zeros(NC.y),S=zeros(NC.x),N=zeros(NC.x)),
        # val     =   (E=400.0*ones(NC.y),W=400.0*ones(NC.y),S=400.0*ones(NC.x),N=400.0*ones(NC.x)),
    )
    # ------------------------------------------------------------------- #
    # Initial Condition ================================================= #
    # Velocity ---
    IniVelocity!(Ini.V,D,VBC,NV,Δ,M,x,y;Ini.ε)
    # Temperature --- 
    IniTemperature!(Ini.T,M,NC,D,x,y;Tb=400.0,Ta=600.0)
    # ------------------------------------------------------------------- #
    # Time ============================================================== #
    T   =   TimeParameter( 
        tmax    =   0.0,  
        Δfacc   =   0.9,                #   Courant time factor
        Δfacd   =   0.9,                #   Diffusion time factor
        itmax   =   200,                #   Maximum iterations
    )
    T.tmax  =   20.589 * 1e6 * (60*60*24*365.25)   # [ s ]
    T.Δc    =   T.Δfacc * minimum((Δ.x,Δ.y)) / 
                    (sqrt(maximum(abs.(D.vx))^2 + maximum(abs.(D.vy))^2))
    T.Δd    =   T.Δfacd * (1.0 / (2.0 * P.κ *(1.0/Δ.x^2 + 1/Δ.y^2)))
    T.Δ     =   minimum([T.Δd,T.Δc])
    Time    =   zeros(T.itmax)
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
    # Numbering, without ghost nodes! ---
    off    = [  NV.x*NC.y,                          # vx
                NV.x*NC.y + NC.x*NV.y,              # vy
                NV.x*NC.y + NC.x*NV.y + NC.x*NC.y]  # Pt

    Num    =    (
        Vx  =   reshape(1:NV.x*NC.y, NV.x, NC.y), 
        Vy  =   reshape(off[1]+1:off[1]+NC.x*NV.y, NC.x, NV.y), 
        Pt  =   reshape(off[2]+1:off[2]+NC.x*NC.y,NC...),#
        T   =   reshape(1:NC.x*NC.y, NC.x, NC.y),
    )
    ndof    =   maximum(Num.T)   
    # Momentum conservation (MCE) ---
    # Iterations --- 
    niterM      =   10
    ϵM          =   1e-8
    KM          =   ExtendableSparseMatrix(ndof,ndof)
    # Energy conservation (ECE) ---
    if FD.Method.Diff==:implicit || FD.Method.Diff==:CN
        if FD.Method.Diff==:CN
            K1      =   ExtendableSparseMatrix(ndof,ndof)
            K2      =   ExtendableSparseMatrix(ndof,ndof)
        else
            K       =   ExtendableSparseMatrix(ndof,ndof)
        end
        rhs         =   zeros(ndof)
    elseif FD.Method.Diff==:dc
        niter       =   10
        ϵ           =   1e-20
        K           =   ExtendableSparseMatrix(ndof,ndof)
        R           =   zeros(Float64,NC...)
        ∂2T         =   (∂x2=zeros(NC.x, NC.y), ∂y2=zeros(NC.x, NC.y),
                        ∂x20=zeros(NC.x, NC.y), ∂y20=zeros(NC.x, NC.y))
    end
    # Residuals ---
    Fm     =    (
        x       =   zeros(Float64,NV.x, NC.y), 
        y       =   zeros(Float64,NC.x, NV.y)
    )
    FPt     =   zeros(Float64,NC...)      
    # ------------------------------------------------------------------- #
    # Time Loop ========================================================= #
    @timeit to "Time Loop" begin
    for it = 1:T.itmax
        δx      =   zeros(maximum(Num.Pt))
        F       =   zeros(maximum(Num.Pt))
        # Update Time ---
        if it>1
            Time[it]  =   Time[it-1] + T.Δ
        end
        @printf("Time step: #%04d, Time [Myr]: %04e\n ",it,
                    Time[it]/(60*60*24*365.25)/1.0e6)
        # Momentum Conservation Equation ================================ #
        # Initial Residual ---
        D.vx[2:end-1,:]    .=  0.0
        D.vy[:,1:end-1]    .=  0.0
        D.Pt               .=  0.0
        @timeit to "Solution Iteration" begin
        @printf("---Momentum Calculation ---\n")        
        for iter = 1:niterM
            @timeit to "Residual" begin
            Residuals2D!(D,VBC,ε,τ,divV,Δ,D.ηc,D.ηv,P.g,Fm,FPt)
            F[Num.Vx]   =   Fm.x[:]
            F[Num.Vy]   =   Fm.y[:]
            F[Num.Pt]   =   FPt[:]
            @printf("||R|| = %1.4e\n", norm(F)/length(F))
            norm(F)/length(F) < ϵM ? break : nothing
            end
            # Assemble Coefficients --
            @timeit to "Assembly" begin
            KM      =   Assembly(NC, NV, Δ, D.ηc, D.ηv, VBC, Num)
            end
            # ---
            # Solution of the linear system --
            @timeit to "Solution" begin
            δx      =   - KM \ F
            end
            # ---
            # Update Unknown Variables ---
            D.vx[:,2:end-1]     .+=  δx[Num.Vx]
            D.vy[2:end-1,:]     .+=  δx[Num.Vy]
            D.Pt                .+=  δx[Num.Pt]
        end
        end
        # Get second invariants ---
        GetSecondInvariants!(ε,τ)
        # Get the velocity on the centroids --- 
        # For visualization purposes only and Euler advection ---
        for i = 1:NC.x
            for j = 1:NC.y
                D.vxc[i,j]  = (D.vx[i,j+1] + D.vx[i+1,j+1])/2
                D.vyc[i,j]  = (D.vy[i+1,j] + D.vy[i+1,j+1])/2
            end
        end
        @. D.vc        = sqrt(D.vxc^2 + D.vyc^2)
        # # ---
        # @show(minimum(D.vc))
        # @show(maximum(D.vc))
        # # ---
        # --------------------------------------------------------------- #
        # Update time stepping ========================================== #
        T.Δc        =   T.Δfacc * minimum((Δ.x,Δ.y)) / 
                (sqrt(maximum(abs.(D.vx))^2 + maximum(abs.(D.vy))^2))
        T.Δd        =   T.Δfacd * (1.0 / (2.0 * P.κ *(1.0/Δ.x^2 + 1/Δ.y^2)))
        T.Δ         =   minimum([T.Δd,T.Δc])
        if Time[it] > T.tmax
            T.Δ         =   T.tmax - Time[it-1]
            Time[it]    =   Time[it-1] + T.Δ
            it          =   T.itmax
        end
        # Energy Conservation Equation ================================== #
        # Shear heating calculation ---
        @. Hₛ   =   τ.xx .* ε.xx .+ τ.yy .* ε.yy .+ 2.0 .* τ.xyc .* ε.xyc
        # ---
        # Temperature diffusion --- 
        @printf("---Energy Calculation---\n")
        if FD.Method.Diff==:explicit
            ForwardEuler2Dc!(D, P.κ, Δ.x, Δ.y, T.Δ, NC, TBC;
                        ρ₀ = P.ρ₀, cp = P.cp, Qₛ = Hₛ)
        elseif FD.Method.Diff==:implicit
            BackwardEuler2Dc!(D, P.κ, Δ.x, Δ.y, T.Δ, NC, TBC, rhs, K, Num; 
                        ρ₀ = P.ρ₀, cp = P.cp, Qₛ = Hₛ)
        elseif FD.Method.Diff==:CN
            CNA2Dc!(D, P.κ, Δ.x, Δ.y, T.Δ, NC, TBC, rhs, K1, K2, Num; 
                        ρ₀ = P.ρ₀, cp = P.cp, Qₛ = Hₛ)
        elseif FD.Method.Diff==:dc
            D.T0    .=  D.T
            for iter = 1:niter
                # Evaluate residual
                ComputeResiduals2Dc!(R, D.T, D.T_ex, D.T0, D.T_ex0, ∂2T, 
                        P.κ, TBC, Δ, T.Δ;
                        C = C, ρ₀ = P.ρ₀, cp = P.cp, Qₛ = Hₛ)
                @printf("||R|| = %1.4e\n", norm(R)/length(R))
                norm(R)/length(R) < ϵ ? break : nothing
                # Assemble linear system
                K  = AssembleMatrix2Dc(P.κ, TBC, Num, NC, Δ, T.Δ;C=C)
                # Solve for temperature correction: Cholesky factorisation
                Kc = cholesky(K.cscmatrix)
                # Solve for temperature correction: Back substitutions
                δT = -(Kc\R[:])
                # Update temperature
                @. D.T += δT[Num.T]
            end
        end
        if FD.Method.Adv==:tracers 
            # Update tracer info from grid ---
        else
            # Advect temperature, Eulerian schemes ---
            # Correct velocity for eulerian advection schemes ---
            @. vx_bg_c  =   -Ini.ε .* x.c2d
            @. vy_bg_c  =   Ini.ε .* y.c2d
            @. vxc_res  =   D.vxc - vx_bg_c 
            @. vyc_res  =   D.vyc - vy_bg_c
            if it == 1
                @. vxco_res    =   D.vxc - vx_bg_c 
                @. vyco_res    =   D.vyc - vy_bg_c 
            end
            if FD.Method.Adv==:upwind
                upwindc2D!(D.T,D.T_ex,vxc_res,vyc_res,NC,T.Δ,Δ.x,Δ.y)
            elseif FD.Method.Adv==:slf
                slfc2D!(D.T,D.T_ex,D.T_ex0,vxc_res,vyc_res,NC,T.Δ,Δ.x,Δ.y)
            elseif FD.Method.Adv==:semilag
                semilagc2D!(D.T,D.T_ex,vxc_res,vyc_res,vxco_res,vyco_res,x,y,T.Δ)
                vxco_res    .=   vxc_res
                vyco_res    .=   vyc_res
            end
        end
        # Advection of tracers ---
        @printf("Running on %d thread(s)\n", nthreads())  
        AdvectTracer2D(Ma,nmark,D,x,y,T.Δ,Δ,NC,rkw,rkv)
        # ---
        # --------------------------------------------------------------- #
        if mod(it,2) == 0 || it == T.itmax || it == 1
            # # r = heatmap(x.c./1e3,y.c./1e3,D.ρ',color=:inferno,
            # #             xlabel="x[km]",ylabel="y[km]",colorbar=false,
            # #             title="ρ",
            # #             aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),                             ylims=(M.ymin/1e3, M.ymax/1e3),
            # #             layout=(3,1),subplot=1)
            # # p = scatter(Ma.x[1:Pl.mainc:end]./1e3,Ma.y[1:Pl.mainc:end]./1e3,
            # #             ms=1,ma=0.5,mc=Ma.phase[1:Pl.mainc:end],markerstrokewidth=0.0,
            # #             xlabel="x[km]",ylabel="y[km]",colorbar=true,
            # #             title="tracers",label="",
            # #             aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3), 
            # #             ylims=(M.ymin/1e3, M.ymax/1e3),
            # #             layout=(2,2),subplot=1)
            # p = heatmap(x.c./1e3,y.c./1e3,D.p',color=:inferno,
            #             xlabel="x[km]",ylabel="y[km]",colorbar=true,
            #             title="phase",
            #             aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),                             
            #             ylims=(M.ymin/1e3, M.ymax/1e3),
            #             layout=(2,2),subplot=1)
            # heatmap!(p,x.c./1e3,y.c./1e3,(D.vc.*Pl.qsc)',
            #             xlabel="x[km]",ylabel="y[km]",colorbar=true,
            #             title="V_c [cm/a]",color=cgrad(:tokyo),
            #             aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),
            #             ylims=(M.ymin/1e3, M.ymax/1e3),
            #             layout=(2,2),subplot=3)
            # # heatmap!(p,x.v./1e3,y.ce./1e3,(D.vx.*Pl.qsc)',
            # #             xlabel="x[km]",ylabel="y[km]",colorbar=true,
            # #             title="V_x [cm/a]",color=cgrad(:batlow),
            # #             aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),
            # #             ylims=(M.ymin/1e3, M.ymax/1e3),
            # #             layout=(2,2),subplot=3)
            # # heatmap!(p,x.ce./1e3,y.v./1e3,(D.vy.*Pl.qsc)',
            # #             xlabel="x[km]",ylabel="y[km]",colorbar=true,
            # #             title="V_y [cm/a]",color=cgrad(:batlow),
            # #             aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),
            # #             ylims=(M.ymin/1e3, M.ymax/1e3),
            # #             layout=(2,2),subplot=4)
            r = heatmap(x.c./1e3,y.c./1e3,log10.(D.ηc'),color=reverse(cgrad(:roma)),
                        xlabel="x[km]",ylabel="y[km]",title="η_c",
                        # clims=(15,27),
                        aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3), 
                        ylims=(M.ymin/1e3, M.ymax/1e3),colorbar=true,
                        layout=(2,2),subplot=1)
            heatmap!(r,x.c./1e3,y.c./1e3,log10.(ε.II'),color=cgrad(:batlow),
                        xlabel="x[km]",ylabel="y[km]",title="ε_II",
                        # clims=(15,27),
                        aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3), 
                        ylims=(M.ymin/1e3, M.ymax/1e3),colorbar=true,
                        layout=(2,2),subplot=3)
            # quiver!(p,x.c2d[1:Pl.qinc:end,1:Pl.qinc:end]./1e3,
            #             y.c2d[1:Pl.qinc:end,1:Pl.qinc:end]./1e3,
            #             quiver=(D.vxc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc,
            #                     D.vyc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc),        
            #             la=0.5,color="white",
            #             layout=(2,2),subplot=3)
            heatmap!(r,x.c./1e3,y.c./1e3,Hₛ',color=:inferno,
                        xlabel="x[km]",ylabel="y[km]",colorbar=true,
                        title="Hₛ",
                        aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),
                        ylims=(M.ymin/1e3, M.ymax/1e3),
                        layout=(2,2),subplot=2)
            heatmap!(r,x.c./1e3,y.c./1e3,D.T',color=:lajolla,
                        xlabel="x[km]",ylabel="y[km]",colorbar=true,
                        title="T",
                        aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),
                        ylims=(M.ymin/1e3, M.ymax/1e3),
                        layout=(2,2),subplot=4)
            if save_fig == 1
                Plots.frame(anim)
            elseif save_fig == 0
                # display(p)
                # display(q)
                display(r)
                # display(t)
                # display(s)
                # display(u)
            end
        end
        # Deform and remesh grid ======================================== #
        # Calculate new xmin, xmax and ymax ---
        M.xmin  =   M.xmin * exp(-Ini.ε * T.Δ)
        M.xmax  =   M.xmax * exp(-Ini.ε * T.Δ)
        M.ymax  =   M.ymax * exp(Ini.ε * T.Δ)
        # Remeshing ---
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
        # --------------------------------------------------------------- #
        # Update boundary velocity ====================================== # 
        IniVelocity!(Ini.V,D,VBC,NV,Δ,M,x,y;Ini.ε)
        # --------------------------------------------------------------- #
        # Update Rheology =============================================== #
        # if FD.Method.Adv==:tracer
            # Update everything
        # else
            # Update only phase 
        # end
        CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,NV)
        @timeit to "Tracer Interpolation" begin
        # Interpolate phase from tracers to grid ---
        Markers2Cells(Ma,nmark,MAVG.PC_th,D.ρ_ex,MAVG.wte_th,D.wte,x,y,Δ,Aparam,ρ)
        D.ρ     .=   D.ρ_ex[2:end-1,2:end-1]  
        Markers2Cells(Ma,nmark,MAVG.PC_th,D.p_ex,MAVG.wte_th,D.wte,x,y,Δ,Aparam,phase)
        D.p     .=  D.p_ex[2:end-1,2:end-1]
        Markers2Cells(Ma,nmark,MAVG.PC_th,D.η_ex,MAVG.wte_th,D.wte,x,y,Δ,Aparam,η)
        D.ηc    .=   D.η_ex[2:end-1,2:end-1]
        Markers2Vertices(Ma,nmark,MAVG.PV_th,D.ηv,MAVG.wtv_th,D.wtv,x,y,Δ,Aparam,η)
        end
        # end
        # Calculate overall strain of the box --- 
        strain      +=   Ini.ε * T.Δ
        @show strain
        if strain >= 0.25
            @printf("Strain %g reached maximum strain of 0.25\n",strain)
            break
        end
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
using Plots, ExtendableSparse
using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
using GeoModBox.AdvectionEquation.TwoD, GeoModBox.HeatEquation.TwoD
# using GeoModBox.Tracers.TwoD
# using Base.Threads
using Statistics, LinearAlgebra
using Printf

function ThermalConvection()
    # Define numerical methods ========================================== #
    # Advection Scheme ---
    #   1) upwind, 2) slf, 3) semilag, 4) tracers
    #       --- slf not working ---
    # Diffusion Scheme --- 
    #   1) explicit, 2) implicit, 3) CNA, 4) ADI, 5) dc
    #       dc - source term missing!
    # Momentum Equation --- 
    #   1) direct, 2) dc 
    FD          =   (Method     = (
        Adv=:semilag,
        Diff=:CNA,
        Mom=:direct),
    )
    # Define Initial Condition ---
    # Temperature - 
    #   1) circle, 2) gaussian, 3) block, 4) linear, 5) lineara
    # !!! Gaussian is not working!!! 
    # Velocity - 
    #   1) RigidBody, 2) ShearCell
    Ini         =   (T=:lineara,) 
    # ------------------------------------------------------------------- #
    # Plot Einstellungen ================================================ #
    Pl  =   (
        qinc        =   5,
        qsc         =   2e2*(100*(60*60*24*365.15)),
    )
    # Animationssettings ================================================ #
    k           =   scatter()
    path        =   string("./exercises/Correction/Results/")
    anim        =   Plots.Animation(path, String[] )
    save_fig    =   1
    # ------------------------------------------------------------------- #
    # Define model geometry constants =================================== #
    M   =   (
        xmin    =   0.0,            #   [ m ] 
        xmax    =   8700e3,         #   [ m ]
        ymin    =   -2900e3,        #   [ m ]
        ymax    =   0.0,    
    )
    # ------------------------------------------------------------------- #
    # Grid ============================================================== #
    NC  =   (
        x   =   150,    #   horizontal grid resolution
        y   =   50,     #   vertical grid resolution
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
    # Define physical constants ========================================= #
    P   =   (
        g   =   9.81,                   #   Schwerebeschleunigung [m/s^2]
        ρ₀  =   3300.0,                 #   Hintergunddichte [kg/m^3]
        k   =   4.125,                  #   Thermische Leitfaehigkeit [ W/m/K ]
        cp  =   1250.0,                 #   Heat capacity [ J/kg/K ]
        α   =   2.0e-5,                 #   Thermischer Expnasionskoef. [ K^-1 ]
        Q₀  =   [0.0],                  #   Waermeproduktionsrate pro Volumen [W/m^3]
        η₀  =   [3.947725485e23],       #   Viskositaet [ Pa*s ] [1.778087025e21]
    )
    P1  =   (
        κ       =   P.k/P.ρ₀/P.cp,      # 	Thermische Diffusivitaet [ m^2/s ]
        # Q₀      =   P.Q₀[1]/P.ρ₀,       #   Waermeproduktionsrate pro Masse [W/kg]
        ΔT      =   2500.0,             #   Temperaturdifferenz
        # Falls Ra < 0 gesetzt ist, dann wird Ra aus den obigen Parametern
        # berechnet. Falls Ra gegeben ist, dann wird die Referenzviskositaet so
        # angepasst, dass die Skalierungsparameter die gegebene Rayleigh-Zahl
        # ergeben.
        Ra      =   [1.0e5],       #   Rayleigh number
    )
    P   =   merge(P,P1)
    filename    =   string("11_ThermalConvection_",P.Ra[1],
                        "_",Ini.T,"_",FD.Method.Adv,"_",FD.Method.Diff,
                        "_",FD.Method.Mom)
    # ------------------------------------------------------------------- #
    # Allocation ======================================================== #
    D       =   (
        Q       =   zeros(Float64,(NC...)),
        T       =   zeros(Float64,(NC...)),
        T0      =   zeros(Float64,(NC...)),
        T_ex    =   zeros(Float64,(NC.x+2,NC.y+2)),
        T_exo   =   zeros(Float64,(NC.x+2,NC.y+2)),
        ρ       =   zeros(Float64,(NC...)),
        cp      =   zeros(Float64,(NC...)),
        vx      =   zeros(Float64,(NV.x,NV.y+1)),
        vy      =   zeros(Float64,(NV.x+1,NV.y)),    
        Pt      =   zeros(Float64,(NC...)),
        vxc     =   zeros(Float64,(NC...)),
        vyc     =   zeros(Float64,(NC...)),
        vc      =   zeros(Float64,(NC...)),
        wt      =   zeros(Float64,(NC...)),
        wtv     =   zeros(Float64,(NV...)),
        ΔTtop   =   zeros(Float64,NC.x),
        ΔTbot   =   zeros(Float64,NC.x),
        Tmax    =   [0.0],
        Tmin    =   [0.0],
        Tmean   =   [0.0],
    )
    tmp     =   zeros(Float64,(NC...))
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
    # Residuals ---
    Fm     =    (
        x       =   zeros(Float64,NV.x, NC.y), 
        y       =   zeros(Float64,NC.x, NV.y)
    )
    FPt         =   zeros(Float64,NC...)
    # ------------------------------------------------------------------- #
    # Anfangsbedingungen ================================================ #
    # Temperature ------
    Ttop    =   273.15
    Tbot    =   Ttop + P.ΔT
    if FD.Method.Adv==:tracers 
        # Tracer Initialization ---
        # Need to implement incremental marker update first! 
        nmx,nmy     =   3,3
        noise       =   0
        nmark       =   nmx*nmy*NC.x*NC.y
        Aparam      =   :thermal
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
        Ma      =   IniTracer2D(Aparam,nmx,nmy,Δ,M,NC,noise,Ini.p,phase)
        # RK4 weights ---
        rkw     =   1.0/6.0*[1.0 2.0 2.0 1.0]   # for averaging
        rkv     =   1.0/2.0*[1.0 1.0 2.0 2.0]   # for time stepping
        # Count marker per cell ---
        CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,NV,1)
        # Interpolate from markers to cell ---
        Markers2Cells(Ma,nmark,MPC.PG_th,D.ρ,MPC.wt_th,D.wt,x,y,Δ,Aparam,ρ)
    else
        IniTemperature!(Ini.T,M,NC,Δ,D,x,y;Tb=Tbot,Ta=Ttop)
        if FD.Method.Adv==:slf
            D.T_exo    .=  D.T_ex
        end
    end
    # Heat production rate ------
    @. D.Q      = P.Q₀
    # Density ------
    # Since we have a Boussinesq approximation, density is the reference 
    # density
    @. D.ρ      =   P.ρ₀
    # ------------------------------------------------------------------- #
    # Boundary Conditions =============================================== #
    # Temperature ------
    TBC     = (
        type    = (W=:Neumann, E=:Neumann,N=:Dirichlet,S=:Dirichlet),
        val     = (W=zeros(NC.y),E=zeros(NC.y),
                        N=273.15.*ones(NC.x),S=P.ΔT.*ones(NC.x)))
    # Velocity ------
    VBC     =   (
        type    =   (E=:freeslip,W=:freeslip,S=:freeslip,N=:freeslip),
        val     =   (E=zeros(NV.y),W=zeros(NV.y),S=zeros(NV.x),N=zeros(NV.x)),
    )
    # ------------------------------------------------------------------- # 
    # Time ============================================================== #
    T   =   (
        year    =   365.25*3600*24,     #   Seconds per year
        tmax    =   [100000.0],         #   [ Ma ]
        Δfacc   =   1.0,                #   Courant time factor
        Δfacd   =   0.9,                #   Diffusion time factor
        Δ       =   [0.0],
        Δc      =   [0.0],              #   Courant time step
        Δd      =   [0.0],              #   Diffusion time stability criterion
        itmax   =   20000,              #   Maximum iterations
    )
    T.tmax[1]   =   T.tmax[1]*1e6*T.year    #   [ s ]
    T.Δc[1]     =   T.Δfacc * minimum((Δ.x,Δ.y)) / 
                        (sqrt(maximum(abs.(D.vx))^2 + maximum(abs.(D.vy))^2))
    T.Δd[1]     =   T.Δfacd * (1.0 / (2.0 * P.κ *(1.0/Δ.x^2 + 1/Δ.y^2)))

    T.Δ[1]      =   minimum([T.Δd[1],T.Δc[1]])
    Time        =   zeros(T.itmax)
    Nus         =   zeros(T.itmax)
    meanV       =   zeros(T.itmax)
    meanT       =   zeros(T.itmax,NC.y+2)
    find        =   0
    # ------------------------------------------------------------------- #
    # Rayleigh number conditions ======================================== #
    if P.Ra[1] < 0
        # If the Rayleigh number is not explicitly defined, it will be 
        # calculated. 
        P.Ra[1]     =   P.ρ₀*P.g*P.α*P.ΔT*(M.ymax-M.ymin)^3/P.η₀[1]/P.κ
    else
        # If the Rayleigh number is defined, we need to adjust a variable, 
        # e.g., the reference viscosity η₀
        P.η₀[1]     =   P.ρ₀*P.g*P.α*P.ΔT*(M.ymax-M.ymin)^3/P.Ra[1]/P.κ
    end
    # =================================================================== #
    # Linear Equations ================================================== #
    # Momentum Equation ======
    off    = [  NV.x*NC.y,                          # vx
                NV.x*NC.y + NC.x*NV.y,              # vy
                NV.x*NC.y + NC.x*NV.y + NC.x*NC.y]  # Pt

    Num    =    (
        Vx  =   reshape(1:NV.x*NC.y, NV.x, NC.y), 
        Vy  =   reshape(off[1]+1:off[1]+NC.x*NV.y, NC.x, NV.y), 
        Pt  =   reshape(off[2]+1:off[2]+NC.x*NC.y,NC...),
        T   =   reshape(1:NC.x*NC.y, NC.x, NC.y),
    )
    ndof    =   maximum(Num.T)        
    # Temperature Equation ======
    if FD.Method.Diff==:implicit || FD.Method.Diff==:CNA
        if FD.Method.Diff==:CNA
            K1      =   ExtendableSparseMatrix(ndof,ndof)
            K2      =   ExtendableSparseMatrix(ndof,ndof)
        else
            K       =   ExtendableSparseMatrix(ndof,ndof)
        end
        rhs         =   zeros(ndof)
    elseif FD.Method.Diff==:dc
        niter       =   10
        ϵ           =   1e-10
        @. D.cp     =   P.cp
        k           =   (x=zeros(NC.x+1,NC.y), y=zeros(NC.x,NC.y+1))
        @. k.x      =   P.k
        @. k.y      =   P.k
        K           =   ExtendableSparseMatrix(ndof,ndof)
        R           =   zeros(Float64,NC...)
        ∂T          =   (∂x=zeros(NC.x+1, NC.y), ∂y=zeros(NC.x, NC.y+1))
        q           =   (x=zeros(NC.x+1, NC.y), y=zeros(NC.x, NC.y+1))
    end
    # ------------------------------------------------------------------- #
    # Time Loop ========================================================= #
    for it = 1:T.itmax
        χ       =   zeros(maximum(Num.Pt))      #   Unknown Vector ME
        rhsM    =   zeros(maximum(Num.Pt))      #   Right-hand Side ME
        if it>1
            Time[it]  =   Time[it-1] + T.Δ[1]
        end
        @printf("Time step: #%04d, Time [Myr]: %04e\n",it,
                        Time[it]/(60*60*24*365.25)/1.0e6)
        # # Momentum Equation(ME) =======
        D.vx    .=  0.0
        D.vy    .=  0.0 
        D.Pt    .=  0.0
        if FD.Method.Mom==:direct
            # Update K ---
            KM      =   Assemblyc(NC, NV, Δ, P.η₀[1], VBC, Num)
            # Update RHS ---
            # rhs term defined by the Boussinesq approximation
            rhsM    =   updaterhsc( NC, NV, Δ, P.η₀[1], -P.ρ₀*P.α.*D.T, -P.g, VBC, Num )
            # Solve System of Equations ---
            χ       =   KM \ rhsM
            # Update Unknown Variables ---
            D.vx[:,2:end-1]     .=  χ[Num.Vx]
            D.vy[2:end-1,:]     .=  χ[Num.Vy]
            D.Pt                .=  χ[Num.Pt]
        elseif FD.Method.Mom==:dc
            # Initial Residual -------------------------------------------------- #
            @. D.ρ  =   -P.ρ₀*P.α*D.T
            Residuals2Dc!(D,VBC,ε,τ,divV,Δ,P.η₀[1],P.g,Fm,FPt)
            rhsM[Num.Vx]    =   Fm.x[:]
            rhsM[Num.Vy]    =   Fm.y[:]
            rhsM[Num.Pt]    =   FPt[:]
            # ------------------------------------------------------------------- #
            # Assemble Coefficients ============================================= #
            K       =   Assemblyc(NC, NV, Δ,P. η₀[1], VBC, Num)
            # ------------------------------------------------------------------- #
            # Solution of the linear system ===================================== #
            χ      =   - K \ rhsM
            # ------------------------------------------------------------------- #
            # Update Unknown Variables ========================================== #
            D.vx[:,2:end-1]     .+=  χ[Num.Vx]
            D.vy[2:end-1,:]     .+=  χ[Num.Vy]
            D.Pt                .+=  χ[Num.Pt]
            @. D.ρ  =   P.ρ₀
        end
        # ======
        # Get the velocity on the centroids ------
        for i = 1:NC.x
            for j = 1:NC.y
                D.vxc[i,j]  = (D.vx[i,j+1] + D.vx[i+1,j+1])/2
                D.vyc[i,j]  = (D.vy[i+1,j] + D.vy[i+1,j+1])/2
            end
        end
        @. D.vc        = sqrt(D.vxc^2 + D.vyc^2)
        # ---
        @show(maximum(D.vc))
        @show(minimum(D.Pt))
        @show(maximum(D.Pt))
        # Calculate time stepping ======================================= #
        T.Δc[1]     =   T.Δfacc * minimum((Δ.x,Δ.y)) / 
                (sqrt(maximum(abs.(D.vx))^2 + maximum(abs.(D.vy))^2))
        T.Δd[1]     =   T.Δfacd * (1.0 / (2.0 * P.κ *(1.0/Δ.x^2 + 1/Δ.y^2)))
        T.Δ[1]      =   minimum([T.Δd[1],T.Δc[1]])
        if Time[it] > T.tmax[1] 
            T.Δ[1]      =   T.tmax[1] - Time[it-1]
            Time[it]    =   Time[it-1] + T.Δ[1]
            it          =   T.itmax
        end
        # Plot ========================================================== #
        if mod(it,10) == 0 || it == T.itmax || it == 1
            p = heatmap(x.c./1e3,y.c./1e3,D.T',
                    xlabel="x[km]",ylabel="y[km]",colorbar=true,
                    title="Temperature",color=cgrad(:lajolla),
                    aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),
                    ylims=(M.ymin/1e3, M.ymax/1e3),
                    layout=(2,1),subplot=1)
            quiver!(p,x.c2d[1:Pl.qinc:end,1:Pl.qinc:end]./1e3,
                y.c2d[1:Pl.qinc:end,1:Pl.qinc:end]./1e3,
                quiver=(D.vxc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc,
                        D.vyc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc),
                la=0.5,color="white",
                layout=(2,1),subplot=1)
            heatmap!(p,x.c./1e3,y.c./1e3,D.vc',color=:imola,
                xlabel="x[km]",ylabel="y[km]",colorbar=true,
                title="Velocity",aspect_ratio=:equal,
                xlims=(M.xmin/1e3, M.xmax/1e3), 
                ylims=(M.ymin/1e3, M.ymax/1e3),
                layout=(2,1),subplot=2)
            if save_fig == 1
                Plots.frame(anim)
            elseif save_fig == 0
                display(p)
            end
        end
        # --------------------------------------------------------------- #
        # Advection ===================================================== #
        if FD.Method.Adv==:upwind
            upwindc2D!(D.T,D.T_ex,D.vxc,D.vyc,NC,T.Δ[1],Δ.x,Δ.y)            
        elseif FD.Method.Adv==:slf
            slfc2D!(D.T,D.T_ex,D.T_exo,D.vxc,D.vyc,NC,T.Δ[1],Δ.x,Δ.y)
        elseif FD.Method.Adv==:semilag
            semilagc2D!(D.T,D.T_ex,D.vxc,D.vyc,[],[],x,y,T.Δ[1])
        elseif FD.Method.Adv==:tracers
            # Advect tracers ---
            @printf("Running on %d thread(s)\n", nthreads())  
            AdvectTracer2D(Ma,nmark,D,x,y,T.Δ[1],Δ,NC,rkw,rkv,1)
            CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,NV,it)
            # Interpolate phase from tracers to grid ---
            Markers2Cells(Ma,nmark,MPC.PG_th,D.ρ,MPC.wt_th,D.wt,x,y,Δ,Aparam,ρ)
        end
        # --------------------------------------------------------------- #
        # Diffusion ===================================================== #
        if FD.Method.Diff==:explicit
            ForwardEuler2Dc!(D, P.κ, Δ.x, Δ.y, T.Δ[1], D.ρ, P.cp, NC, TBC)
        elseif FD.Method.Diff==:implicit
            BackwardEuler2Dc!(D, P.κ, Δ.x, Δ.y, T.Δ[1], D.ρ, P.cp, NC, TBC, rhs, K, Num)
        elseif FD.Method.Diff==:CNA
            CNA2Dc!(D, P.κ, Δ.x, Δ.y, T.Δ[1], D.ρ, P.cp, NC, TBC, rhs, K1, K2, Num)
        elseif FD.Method.Diff==:ADI
            ADI2Dc!(D, P.κ, Δ.x, Δ.y, T.Δ[1], D.ρ, P.cp, NC, TBC)
        elseif FD.Method.Diff==:dc
            D.T0    .=  D.T
            for iter = 1:niter
                # Evaluate residual
                ComputeResiduals2D!(R, D.T, D.T_ex, D.T0, ∂T, q, D.ρ, D.cp, k, TBC, Δ, T.Δ[1])
                # @printf("||R|| = %1.4e\n", norm(R)/length(R))
                norm(R)/length(R) < ϵ ? break : nothing
                # Assemble linear system
                K  = AssembleMatrix2D(D.ρ, D.cp, k, TBC, Num, NC, Δ, T.Δ[1])
                # Solve for temperature correction: Cholesky factorisation
                Kc = cholesky(K.cscmatrix)
                # Solve for temperature correction: Back substitutions
                δT = -(Kc\R[:])
                # Update temperature
                @. D.T += δT[Num.T]
            end        
        end
        # --------------------------------------------------------------- #
        # Heat flow at the surface ====================================== #
        @. D.ΔTbot  =   
                (((D.T_ex[2:end-1,2]+D.T_ex[2:end-1,3])/2.0) - 
                ((D.T_ex[2:end-1,2]+D.T_ex[2:end-1,1])/2.0)) / Δ.y
        @. D.ΔTtop  =   
                (((D.T_ex[2:end-1,end-2]+D.T_ex[2:end-1,end-1]) / 2.0) - 
                ((D.T_ex[2:end-1,end-1]+D.T_ex[2:end-1,end]) / 2.0)) / Δ.y
        Nus[it]     =   mean(D.ΔTtop)
        meanT[it,:] =   mean(D.T_ex,dims=1)
        meanV[it]   =   mean(D.vc)
        # --------------------------------------------------------------- #
        # Check break =================================================== #
        # If the maximum time is reached or if the models reaches steady, 
        # state the time loop is stoped! 
        if Time[it]/1e6/T.year > 50     # [ Ma ]
            epsC    =   1e-16; 
            ind     =   findfirst(Time./1e6/T.year .> 
                            (Time[it]/1e6/T.year - 50))
            epsV    =   std(meanV[ind:it])
            # @show epsV1, epsV
            if save_fig == 1
                @show it,log10((epsV))
                plot!(k,(Time[it]/1e6/T.year,log10((epsV))),
                    xlabel="it",ylabel="log₁₀(εᵥ)",label="",
                    markershape=:circle,markercolor=:black)
            end
            find    =   it
            @printf("ε_V = %g, ε_C = %g \n",epsV,epsC)
            if Time[it] >= T.tmax[1]
                @printf("Maximum time reached!\n")
                find    =   it
                break
            elseif (epsV <= epsC)
                @printf("Convection reaches steady state!\n")
                find    =   it
                break
            end
        end
        # --------------------------------------------------------------- #
        @printf("\n")
    end
    if save_fig == 1
        # Write the frames to a GIF file
        Plots.gif(anim, string( path, filename, ".gif" ), fps = 15)
        foreach(rm, filter(startswith(string(path,"00")), readdir(path,join=true)))
    end
    # Save final figure ===================================================== #
    p2 = heatmap(x.c./1e3,y.c./1e3,D.T',
            xlabel="x[km]",ylabel="y[km]",colorbar=true,
            title="Temperature",color=cgrad(:lajolla),
            aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),
            ylims=(M.ymin/1e3, M.ymax/1e3),
            layout=(2,1),subplot=1)
        quiver!(p2,x.c2d[1:Pl.qinc:end,1:Pl.qinc:end]./1e3,
                y.c2d[1:Pl.qinc:end,1:Pl.qinc:end]./1e3,
                quiver=(D.vxc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc,
                        D.vyc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc),
                la=0.5,color="white",
                layout=(2,1),subplot=1)
        heatmap!(p2,x.c./1e3,y.c./1e3,D.vc',color=:imola,
                xlabel="x[km]",ylabel="y[km]",colorbar=true,
                title="Velocity",aspect_ratio=:equal,
                xlims=(M.xmin/1e3, M.xmax/1e3), 
                ylims=(M.ymin/1e3, M.ymax/1e3),
                layout=(2,1),subplot=2)
    if save_fig == 1
        savefig(k,string("./exercises/Correction/Results/11_ThermalConvection_iterations_",P.Ra[1],
            "_",Ini.T,"_",FD.Method.Adv,"_",FD.Method.Diff,"_",FD.Method.Mom,".png"))
        savefig(p2,string("./exercises/Correction/Results/11_ThermalConvection_",P.Ra[1],
            "_Final_Stage_",Ini.T,"_",FD.Method.Adv,"_",FD.Method.Diff,"_",FD.Method.Mom,".png"))
    else
        display(p2)
    end
    # ----------------------------------------------------------------------- #
    # Plot time serieses ==================================================== #
    q2  =   plot(Time[1:find]./1e6/T.year,Nus[1:find],
                    xlabel="Time [ Ma ]", ylabel="Nus",label="",
                    layout=(2,1),suplot=1)
    plot!(q2,Time[1:find]./1e6/T.year,meanV[1:find],
                xlabel="Time [ Ma ]", ylabel="V_{RMS}",label="",
                layout=(2,1),subplot=2)
    if save_fig == 1
        savefig(q2,string("./exercises/Correction/Results/11_TimeSeries_",P.Ra[1],"_",
                            Ini.T,"_",FD.Method.Adv,"_",FD.Method.Diff,"_",FD.Method.Mom,".png"))
    else
        display(q2)
    end
# ======================================================================= #
end

ThermalConvection()

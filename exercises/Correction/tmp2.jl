using Plots, ExtendableSparse
using GeoModBox
using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
using GeoModBox.AdvectionEquation.TwoD, GeoModBox.HeatEquation.TwoD
using GeoModBox.Scaling
using GeoModBox.Tracers.TwoD
using Base.Threads
using Statistics, LinearAlgebra
using Printf, LaTeXStrings, Measures

function main()
    # Define numerical methods ========================================== #
    # Diffusion Scheme --- 
    #   1) explicit, 2) implicit, 3) CNA, 4) ADI, 5) dc
    # Advection Scheme ---
    #   1) upwind, 2) semilag, 3) tracers
    #       --- tracers need to be modified ---
    # Momentum Equation --- 
    #   1) direct, 2) dc 
    FD          =   (Method     = (
        Diff=:CNA,
        Adv=:upwind,
        Mom=:dc),
    )
    # Define Initial Condition ---
    # Temperature - 
    Ini         =   (T=:lineara,)
    # ------------------------------------------------------------------- #
    # Plot Einstellungen ================================================ #
    Pl  =   (
        qinc        =   5,
        qsc         =   1.0e-4
    )
    # Animationssettings ================================================ #
    # ks          =   scatter()
    path        =   string("./exercises/Correction/Results/")
    anim        =   Plots.Animation(path, String[] )
    save_fig    =   1
    # ------------------------------------------------------------------- #
    # Modellgeometrie Konstanten ======================================== #
    M   =   Geometry(
        xmin    =   0.0,                #   [ m ] 
        xmax    =   8700e3,             #   [ m ]
        ymin    =   -2900e3,            #   [ m ]
        ymax    =   0.0,                #   [ m ]
    )
    # ------------------------------------------------------------------- #
    # Referenzparameter ================================================= #
    P   =   Physics(
        g       =   9.81,               #   Schwerebeschleunigung [m/s^2]
        ρ₀      =   3300.0,             #   Hintergunddichte [kg/m^3]
        k       =   4.125,              #   Thermische Leitfaehigkeit [ W/m/K ]
        cp      =   1250.0,             #   Heat capacity [ J/kg/K ]
        α       =   2.0e-5,             #   Thermischer Expnasionskoef. [ K^-1 ]
        Q₀      =   0.0,                #   Waermeproduktionsrate pro Volumen [W/m^3]
        η₀      =   3.947725485e23,     #   Viskositaet [ Pa*s ] [1.778087025e21]
        ΔT      =   2500.0,             #   Temperaturdifferenz
        # Falls Ra < 0 gesetzt ist, dann wird Ra aus den obigen Parametern
        # berechnet. Falls Ra gegeben ist, dann wird die Referenzviskositaet so
        # angepasst, dass die Skalierungsparameter die gegebene Rayleigh-Zahl
        # ergeben.
        Ra      =   1.0e6,              #   Rayleigh number
        Ttop    =   273.15,             #   Temperatur an der Oberfläche [ K ]
    )
    # ------------------------------------------------------------------- #
    # Rayleigh Zahl Bedingungen ========================================= #
    if P.Ra < 0
        # Falls die Rayleigh Zahl nicht explizit angegeben wird, dann 
        # wird sie hier berechnet
        P.Ra     =   P.ρ₀*P.g*P.α*P.ΔT*(M.ymax-M.ymin)^3/P.η₀/P.κ
    else
        # Falls die Rayleigh Zahl explizit angegeben ist, dann wird hier 
        # die Referenzviskositaet η₀ angepasst. 
        P.η₀     =   P.ρ₀*P.g*P.α*P.ΔT*(M.ymax-M.ymin)^3/P.Ra/P.κ
    end
    # =================================================================== #
    # Definiere Skalierungskonstanten =================================== # 
    S   =   ScalingConstants!(M,P)
    # ------------------------------------------------------------------- #
    # Gittereinstellungen =============================================== #
    NC  =   (
        x   =   150,
        y   =   50,
    )
    NV      =   (
        x   =   NC.x + 1,
        y   =   NC.y + 1,
    )
    Δ       =   GridSpacing(
        x   =   (M.xmax - M.xmin)/NC.x,
        y   =   (M.ymax - M.ymin)/NC.y,
    )
    # ------------------------------------------------------------------- #
    # Initialisierung der Datenfelder =================================== #
    D       =   DataFields(
        Q       =   zeros(Float64,(NC...)),
        T       =   zeros(Float64,(NC...)),
        T0      =   zeros(Float64,(NC...)),
        T_ex    =   zeros(Float64,(NC.x+2,NC.y+2)),
        T_ex0   =   zeros(Float64,(NC.x+2,NC.y+2)),
        Told_ex =   zeros(Float64,(NC.x+2,NC.y+2)),
        ρ       =   ones(Float64,(NC...)),
        vx      =   zeros(Float64,(NV.x,NV.y+1)),
        vy      =   zeros(Float64,(NV.x+1,NV.y)),    
        Pt      =   zeros(Float64,(NC...)),
        vxc     =   zeros(Float64,(NC...)),
        vyc     =   zeros(Float64,(NC...)),
        vxco    =   zeros(Float64,(NC...)),
        vyco    =   zeros(Float64,(NC...)),
        vc      =   zeros(Float64,(NC...)),
        wt      =   zeros(Float64,(NC...)),
        wte     =   zeros(Float64,(NC.x+2,NC.y+2)),
        wtv     =   zeros(Float64,(NV...)),
        ΔTtop   =   zeros(Float64,NC.x),
        ΔTbot   =   zeros(Float64,NC.x),
        Tmax    =   0.0,
        Tmin    =   0.0,
        Tmean   =   0.0,
    )
    # Wärmeproduktionsrate ------
    @. D.Q      =   P.Q₀
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
    # Residuals ---
    Fm     =    (
        x       =   zeros(Float64,NV.x, NC.y), 
        y       =   zeros(Float64,NC.x, NV.y)
    )
    FPt         =   zeros(Float64,NC...)
    # ------------------------------------------------------------------- #
    # Heat Flux calculations ---
    Tv1         =   zeros(NV.x,1)
    Tv2         =   zeros(NV.x,1)
    dTdy        =   zeros(NV.x,1)
    # ------------------------------------------------------------------- #
    # Zeitparameter ===================================================== #
    T   =   TimeParameter(
        tmax    =   1000000.0,          #   [ Ma ]
        Δfacc   =   0.9,                #   Courant time factor
        Δfacd   =   0.9,                #   Diffusion time factor
        itmax   =   10000,              #   Maximum iterations
    )
    T.tmax      =   T.tmax*1e6*T.year    #   [ s ]
    T.Δc        =   T.Δfacc * minimum((Δ.x,Δ.y)) / 
                        (sqrt(maximum(abs.(D.vx))^2 + maximum(abs.(D.vy))^2))
    T.Δd        =   T.Δfacd * (1.0 / (2.0 * P.κ *(1.0/Δ.x^2 + 1/Δ.y^2)))

    T.Δ         =   minimum([T.Δd,T.Δc])
    # Statistics -------------------------------------------------------- #
    Time            =   zeros(T.itmax)
    Nus             =   zeros(T.itmax)
    meanV           =   zeros(T.itmax)
    meanT           =   zeros(T.itmax)
    epsV_history    =   fill(NaN, T.itmax)
    find            =   0
    count           =   0
    epsV0           =   0
    final_step      =   0 
    # ------------------------------------------------------------------- #
    # Skalierungsgesetze ================================================ #
    ScaleParameters!(S,M,Δ,T,P,D)
    # ------------------------------------------------------------------- #
    # Koordinaten ======================================================= #
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
    # Anfangsbedingungen ================================================ #
    # Temperatur ------
    IniTemperature!(Ini.T,M,NC,D,x,y;Tb=P.Tbot,Ta=P.Ttop)
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
        MAVG        = (
            PC_th   =   [similar(D.wte) for _ = 1:nthreads()],  # per thread
            PV_th   =   [similar(D.wtv) for _ = 1:nthreads()],   # per thread
            wte_th  =   [similar(D.wte) for _ = 1:nthreads()],  # per thread
            wtv_th  =   [similar(D.wtv) for _ = 1:nthreads()],  # per thread
        )
        Ma      =   IniTracer2D(Aparam,nmx,nmy,Δ,M,NC,noise,:nothing,:nothing)
        # RK4 weights ---
        rkw     =   1.0/6.0*[1.0 2.0 2.0 1.0]   # for averaging
        rkv     =   1.0/2.0*[1.0 1.0 2.0 2.0]   # for time stepping
        # Count marker per cell ---
        CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,NV)
        # Interpolate from markers to cell ---
        @threads for k = 1:nmark
            Ma.T[k] =   FromCtoM(D.T_ex, k, Ma, x, y, Δ, NC)
        end
        ΔT_grid     =   zeros(Float64,(NC.x+2,NC.y+2))
    end
    # ------------------------------------------------------------------- #
    # Boundary Conditions =============================================== #
    # Temperature ------
    TBC     = (
        type    = (W=:Neumann, E=:Neumann,N=:Dirichlet,S=:Dirichlet),
        val     = (W=zeros(NC.y),E=zeros(NC.y),
                        N=P.Ttop.*ones(NC.x),S=P.Tbot.*ones(NC.x)))
    # Velocity ------
    VBC     =   (
        type    =   (E=:freeslip,W=:freeslip,S=:freeslip,N=:freeslip),
        val     =   (E=zeros(NV.y),W=zeros(NV.y),S=zeros(NV.x),N=zeros(NV.x),
                    vxE=zeros(NC.y),vxW=zeros(NC.y),vyS=zeros(NC.x),vyN=zeros(NC.x)),
    )
    # ------------------------------------------------------------------- # 
    filename    =   string("12_ThermalConvection_",P.Ra,
                            "_",NC.x,"_",NC.y,
                            "_",Ini.T,"_",FD.Method.Adv,"_",FD.Method.Diff,
                            "_",FD.Method.Mom)
    # Lineares Gleichungssystem ========================================= #
    # Impulserhaltung (IEG) ------
    off     = [ NV.x*NC.y,                          # vx
                NV.x*NC.y + NC.x*NV.y,              # vy
                NV.x*NC.y + NC.x*NV.y + NC.x*NC.y ] # Pt

    Num    =    (
        Vx  =   reshape(1:NV.x*NC.y, NV.x, NC.y), 
        Vy  =   reshape(off[1]+1:off[1]+NC.x*NV.y, NC.x, NV.y), 
        Pt  =   reshape(off[2]+1:off[2]+NC.x*NC.y,NC...),
        T   =   reshape(1:NC.x*NC.y, NC.x, NC.y),
    )
    χ       =   zeros(maximum(Num.Pt))      #   Unknown Vector ME
    rhsM    =   zeros(maximum(Num.Pt))      #   Right-hand Side ME
    if FD.Method.Mom==:dc
        ndofM       =   maximum(Num.Pt) 
        niterM      =   50
        atolM       =   1e-10       #   Absolute tolerance
        rtolM       =   1e-8        #   Relative tolerance; r = atolM0/atolM
        RM          =   0.0         #   Initialize absolute residual    
        RMrel       =   0.0         #   Initialize relative residual 
        KM          =   ExtendableSparseMatrix(ndofM,ndofM)
    end
    KM              =   Assemblyc(NC, NV, Δ, 1.0, VBC, Num)
    KMfac           =   lu(KM.cscmatrix)
    # === Temperature Equation ======
    ndofT           =   maximum(Num.T)
    if FD.Method.Diff==:implicit || FD.Method.Diff==:CNA
        if FD.Method.Diff==:CNA
            K1      =   ExtendableSparseMatrix(ndofT,ndofT)
            K2      =   ExtendableSparseMatrix(ndofT,ndofT)
        else
            K       =   ExtendableSparseMatrix(ndofT,ndofT)
        end
        rhs         =   zeros(ndofT)
    elseif FD.Method.Diff==:dc
        niterT      =   10
        ϵT          =   1e-10
        K           =   ExtendableSparseMatrix(ndofT,ndofT)
        R           =   zeros(Float64,NC...)
        δT          =   zeros(ndofT)
        ∂2T         =   (∂x2=zeros(NC.x, NC.y), ∂y2=zeros(NC.x, NC.y),
                        ∂x20=zeros(NC.x, NC.y), ∂y20=zeros(NC.x, NC.y))
    end
    # ------------------------------------------------------------------- #
    # Zeitschleife ====================================================== #
    for it = 1:T.itmax
        R0      =   0
        verbose_step    =   mod(it, 250) == 0 || it == 1 || final_step == 1
        if it>1
            Time[it]  =   Time[it-1] + T.Δ
        end
        verbose_step && @printf("Time step: #%04d, Time [non-dim]: %04e\n ",it,
                        Time[it])
        # IEG ------
        verbose_step && @printf("---Momentum Calculation ---\n")
        if it == 1
            D.vx    .=  0.0
            D.vy    .=  0.0 
            D.Pt    .=  0.0
        end
        if FD.Method.Mom==:direct
            # Update K ---
            # KM      =   Assemblyc(NC, NV, Δ, 1.0, VBC, Num)
            # Update RHS ---
            # rechte Seite definiert durch die Boussinesq Annäherung
            rhsM    =   updaterhsc( NC, NV, Δ, 1.0, -P.Ra*D.T, 1.0, VBC, Num )
            # Lösen des linearen Gleichungssystems ---
            χ       =   KMfac \ rhsM
            # Update unbekannte Variablen ---
            D.vx[:,2:end-1]     .=  χ[Num.Vx]
            D.vy[2:end-1,:]     .=  χ[Num.Vy]
            D.Pt                .=  χ[Num.Pt]
        elseif FD.Method.Mom==:dc
            @. D.ρ  =   -P.Ra*D.T
            # Anfangsresiduum ------
            for iter = 1:niterM
                Residuals2Dc!(D,VBC,ε,τ,divV,Δ,1.0,1.0,Fm,FPt)
                rhsM[Num.Vx]    .=   Fm.x
                rhsM[Num.Vy]    .=   Fm.y
                rhsM[Num.Pt]    .=   FPt
                RM          =   norm(rhsM)/length(rhsM)
                if iter == 1
                    R0 = max(RM, eps())
                end
                RMrel       =   RM/R0
                if verbose_step
                    @printf("   MCE %2d: ||R|| = %1.4e, ||R||/||R₀|| = %1.4e\n",iter,RM,RMrel)
                end
                (RM < atolM || RM/R0 < rtolM) && break
                # verbose_step && @printf("||R_M|| = %1.4e\n", norm(rhsM)/sqrt(length(rhsM)))
                # norm(rhsM)/sqrt(length(rhsM)) < ϵM ? break : nothing
                # Update K ------
                # KM      =   Assemblyc(NC, NV, Δ, 1.0, VBC, Num)
                # Lösen des lineare Gleichungssystems ------
                χ      .=   - ( KMfac \ rhsM )
                # Update unbekante Variablen ------
                D.vx[:,2:end-1]     .+=  χ[Num.Vx]
                D.vy[2:end-1,:]     .+=  χ[Num.Vy]
                D.Pt                .+=  χ[Num.Pt]
            end
        end
        # ======
        # Berechnung der Geschwindikeit auf den Centroids ------
        for i = 1:NC.x
            for j = 1:NC.y
                D.vxc[i,j]  = (D.vx[i,j+1] + D.vx[i+1,j+1])/2
                D.vyc[i,j]  = (D.vy[i+1,j] + D.vy[i+1,j+1])/2
            end
        end
        @. D.vc        = sqrt(D.vxc^2 + D.vyc^2)
        # ---
        # Heat flow at the surface ====================================== #
        # Get temperature at the vertices 
        @. Tv1  =   (D.T_ex[1:end-1,end] + D.T_ex[2:end,end] + 
                        D.T_ex[1:end-1,end-1] + D.T_ex[2:end,end-1])/4
        @. Tv2  =   (D.T_ex[1:end-1,end-1] + D.T_ex[2:end,end-1] + 
                        D.T_ex[1:end-1,end-2] + D.T_ex[2:end,end-2])/4
        # Calculate temperature gradient --- 
        @. dTdy =   -(Tv1 - Tv2)/Δ.y
        # Calculate Nusselt number ---
        Nus[it]     =   0.0
        # Trapezoidal integration -
        for i = 1:NV.x
            if i == 1 || i == NV.x
                afac = 1
            else
                afac = 2
            end
            Nus[it]     += afac * dTdy[i]
        end
        Nus[it]     *=   Δ.x/2
        Nus[it]     /= (M.xmax - M.xmin)
        # Mean Temperature ---
        meanT[it]   =   mean(D.T)
        # Root Mean Square Velocity ---
        meanV[it] = sqrt(mean(D.vxc.^2 .+ D.vyc.^2))
        # --------------------------------------------------------------- #
        # Berechnung der Zeitschrittlänge =============================== #
        T.Δc        =   T.Δfacc * minimum((Δ.x,Δ.y)) / 
                (sqrt(maximum(abs.(D.vx))^2 + maximum(abs.(D.vy))^2))
        T.Δd        =   T.Δfacd * (1.0 / (2.0 *(1.0/Δ.x^2 + 1/Δ.y^2)))
        T.Δ         =   minimum([T.Δd,T.Δc])
        if Time[it] + T.Δ >= T.tmax
            T.Δ         =   T.tmax - Time[it]
            final_step  =   1
        end
        # Plot ========================================================== #
        if verbose_step
            p = plot(layout=(2,1),size=(1200,800),dpi = 300,
                guidefontsize = 22, tickfontsize = 22,
                titlefontsize = 22)
            heatmap!(p,x.c,y.c,D.T',
                    xlabel="",ylabel=L"y",colorbar=true,
                    framestyle =:box,
                    xformatter = _ -> "",
                    left_margin = 10mm,
                    right_margin = 10mm, 
                    title=L"Temperature",color=cgrad(:lajolla),
                    aspect_ratio=:equal,xlims=(M.xmin, M.xmax),
                    ylims=(M.ymin, M.ymax),clims=(0,1),
                    subplot=1)
            heatmap!(p,x.c,y.c,D.vc',color=:imola,
                xlabel=L"x",ylabel=L"y",colorbar=true,
                framestyle =:box,
                title=L"Velocity",aspect_ratio=:equal,
                xlims=(M.xmin, M.xmax), 
                ylims=(M.ymin, M.ymax),
                subplot=2)
            contour!(p,x.c,y.c,D.T',lw=1,
                color="white",cbar=true,alpha=0.5,
                subplot=1)     
            quiver!(p,x.c2d[1:Pl.qinc:end,1:Pl.qinc:end],
                y.c2d[1:Pl.qinc:end,1:Pl.qinc:end],
                quiver=(D.vxc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc,
                        D.vyc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc),
                la=0.5,color="black",
                subplot=2)
            if save_fig == 1
                Plots.frame(anim)
            elseif save_fig == 0
                display(p)
            end
        end
        # --------------------------------------------------------------- #
        # Advektion ===================================================== #
        if FD.Method.Adv==:upwind
            upwindc2D!(D.T,D.T_ex,D.vxc,D.vyc,NC,T.Δ,Δ.x,Δ.y)
        elseif FD.Method.Adv==:semilag
            if it == 1
                D.vxco  =   copy(D.vxc)
                D.vyco  =   copy(D.vyc)
            end
            semilagc2D!(D.T,D.T_ex,D.vxc,D.vyc,D.vxco,D.vyco,x,y,T.Δ)
            D.vxco  =   copy(D.vxc)
            D.vyco  =   copy(D.vyc)
        elseif FD.Method.Adv==:tracers
            # Advect tracers ---
            @printf("Running on %d thread(s)\n", nthreads())  
            @. ΔT_grid     =   D.T_ex - D.Told_ex
            @threads for k = 1:nmark
                local ΔTm       =   FromCtoM(ΔT_grid, k, Ma, x, y, Δ, NC)
                Ma.T[k]     += ΔTm
            end
            AdvectTracer2D(Ma,nmark,D,x,y,T.Δ,Δ,NC,rkw,rkv)
            CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,NV)
            # Interpolate phase from tracers to grid ---
            Markers2Cells(Ma,nmark,MAVG.PC_th,D.T_ex,MAVG.wte_th,D.wte,x,y,Δ,Aparam,0)           
            D.T         .=  D.T_ex[2:end-1,2:end-1]
            D.Told_ex   .=  D.T_ex
        end
        # --------------------------------------------------------------- #
        # Diffusion ===================================================== #
        verbose_step && @printf("---Energy Calculation---\n")
        if FD.Method.Diff==:explicit
            ForwardEuler2Dc!(D, 1.0, Δ.x, Δ.y, T.Δ, NC, TBC)
        elseif FD.Method.Diff==:implicit
            BackwardEuler2Dc!(D, 1.0, Δ.x, Δ.y, T.Δ, NC, TBC, rhs, K, Num)
        elseif FD.Method.Diff==:CNA
            CNA2Dc!(D, 1.0, Δ.x, Δ.y, T.Δ, NC, TBC, rhs, K1, K2, Num)
        elseif FD.Method.Diff==:ADI
            ADI2Dc!(D, 1.0, Δ.x, Δ.y, T.Δ, NC, TBC)
        elseif FD.Method.Diff==:dc
            D.T0    .=  D.T
            # Assemble linear system
            K  = AssembleMatrix2Dc(1.0, TBC, Num, NC, Δ, T.Δ;C=0.5)
            # Solve for temperature correction: Cholesky factorisation
            Kc = cholesky(K.cscmatrix)
            for iter = 1:niterT
                # Evaluate residual
                ComputeResiduals2Dc!( R, D.T, D.T_ex, D.T0, D.T_ex0, ∂2T, 
                    1.0, TBC, Δ, T.Δ; C = 0.5 )
                verbose_step && @printf("||R_T|| = %1.4e\n", norm(R)/sqrt(length(R)))
                norm(R)/sqrt(length(R)) < ϵT ? break : nothing
                # Solve for temperature correction: Back substitutions
                δT .= -(Kc\R[:])
                # Update temperature
                @. D.T += δT[Num.T]
            end        
            # Update extended field for advection scheme --- 
            D.T_ex[2:end-1,2:end-1]     .=  D.T
        end
        # D.Told_ex .= D.T_ex
        # --------------------------------------------------------------- #
        # Check break =================================================== #
        # If the maximum time is reached or if the models reaches steady 
        # state the time loop is stoped! 
        if Time[it] > 0.0038
        epsC    =   0.0001
            ind = searchsortedfirst(@view(Time[1:it]), Time[it] - 0.0038)
            epsV    =   std(@view meanV[ind:it])
            if count == 0
                epsV0   = max(epsV, eps(Float64))
                count   += 1
            end
            epsV        /= epsV0
            epsV_history[it] = epsV
            find    =   it
            verbose_step && @printf("ε_Vr = %g, ε_Cr = %g \n",epsV,epsC)
            if Time[it] >= T.tmax
                @printf("Maximum time reached!\n")
                find    =   it
                break
            elseif (epsV <= epsC)
                @printf("Convection reaches steady state at %i!\n",it)
                find    =   it
                break
            end
        end
        # --------------------------------------------------------------- #
        verbose_step && @printf("\n")
    end
    if save_fig == 1
        # Write the frames to a GIF file
        Plots.gif(anim, string( path, filename, ".gif" ), fps = 15)
        foreach(rm, filter(startswith(string(path,"00")), readdir(path,join=true)))
    end
    valid   =   findall(isfinite, epsV_history)
    ks  = plot(
            valid,
            log10.(epsV_history[valid]),
            xlabel = "it",
            ylabel = "log₁₀(εᵥ/εᵥ₀)",
            label = "",
            markershape = :circle,
            markercolor = :black,
    )
    # Save final figure ===================================================== #
    p2 = plot(layout=(2,1),size=(1200,800),dpi = 300,
                guidefontsize = 22, tickfontsize = 22,
                titlefontsize = 22)
    heatmap!(p2,x.c,y.c,D.T',
            xlabel="",ylabel=L"y",colorbar=true,
            framestyle =:box,
            xformatter = _ -> "",
            left_margin = 10mm,
            right_margin = 10mm, 
            title=L"Temperature",color=cgrad(:lajolla),
            aspect_ratio=:equal,xlims=(M.xmin, M.xmax),
            ylims=(M.ymin, M.ymax),clims=(0,1),
            subplot=1)
    heatmap!(p2,x.c,y.c,D.vc',color=:imola,
        xlabel=L"x",ylabel=L"y",colorbar=true,
        framestyle =:box,
        title=L"Velocity",aspect_ratio=:equal,
        xlims=(M.xmin, M.xmax), 
        ylims=(M.ymin, M.ymax),
        subplot=2)
    contour!(p2,x.c,y.c,D.T',lw=1,
        color="white",cbar=true,alpha=0.5,
        subplot=1)     
    quiver!(p2,x.c2d[1:Pl.qinc:end,1:Pl.qinc:end],
        y.c2d[1:Pl.qinc:end,1:Pl.qinc:end],
        quiver=(D.vxc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc,
                D.vyc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc),
        la=0.5,color="black",
        subplot=2)
    if save_fig == 1
        savefig(ks,string("./exercises/Correction/Results/12_ThermalConvection_Scaled_iterations_",P.Ra,
                "_",NC.x,"_",NC.y,
                "_",Ini.T,"_",FD.Method.Adv,"_",FD.Method.Diff,"_",FD.Method.Mom,".png"))
        savefig(p2,string("./exercises/Correction/Results/12_ThermalConvection_Scaled_Final_Stage_",P.Ra,
                "_",NC.x,"_",NC.y,"_it_",find,"_",
                Ini.T,"_",FD.Method.Adv,"_",FD.Method.Diff,"_",FD.Method.Mom,".png"))
    elseif save_fig == 0
        display(p2)
        display(ks)
    end
    # ----------------------------------------------------------------------- #
    # Plot time serieses ==================================================== #
    q2  =   plot(layout=(2,1),size=(1200,800),dpi = 300,
                    guidefontsize = 22, tickfontsize = 22,
                    titlefontsize = 22,
                    left_margin = 10mm, right_margin = 10mm, 
                    bottom_margin = 10mm, top_margin = 10mm)
    plot!(q2,Time[1:find],Nus[1:find],
                xlabel="", ylabel=L"Nus",label="",
                xformatter = _ -> "",
                framestyle =:box,linewidth=2,
                xlims=(0.0,Time[find]),
                subplot=1)
    plot!(q2,Time[1:find],meanV[1:find],
                xlabel=L"Time [ non-dim ]", ylabel=L"V_{RMS}",
                label="",
                framestyle =:box,linewidth=2,
                xlims=(0.0,Time[find]),
                subplot=2)
    if save_fig == 1
        savefig(q2,string("./exercises/Correction/Results/12_ThermalConvectionTimeSeries_Scaled_",P.Ra,
                            "_",NC.x,"_",NC.y,"_",Ini.T,"_",FD.Method.Adv,"_",FD.Method.Diff,"_",FD.Method.Mom,".png"))
    elseif save_fig == 0
        display(q2)
    end
    # ======================================================================= #
end



start=time()
main()
stop=time()
@printf("Runtime: %g s",stop-start)
using Plots, Interpolations, ExtendableSparse, LinearAlgebra
using GeoModBox.AdvectionEquation.TwoD, GeoModBox.InitialCondition
using GeoModBox.HeatEquation.TwoD

function EnergyEquation()

    # Definition numerischer Verfahren ================================== #
    # Define Advection Scheme ---
    #   1) upwind, 2) slf, 3) semilag
    # Define Diffusion Scheme --- 
    #   1) explicit, 2) implicit, 3) CNA, 4) ADI, 5) dc
    FD          =   (Method     = (Adv=:semilag,Diff=:none),)
    # Define Initial Condition ---
    # Temperature - 
    #   1) circle, 2) gaussian, 3) block
    # Velocity - 
    #   1) RigidBody, 2) ShearCell
    Ini         =   (T=:circle,V=:RigidBody,) 
    # ------------------------------------------------------------------- #
    # Plot Einstellungen ================================================ #
    Pl  =   (
        inc         =   5,
        sc          =   20*(100*(60*60*24*365.15)),
    )
    # ------------------------------------------------------------------- #
    # Model Konstanten ================================================== #
    M   =   (
        xmin    =   0.0,
        xmax    =   200.0e3,
        ymin    =   0.0,
        ymax    =   200.0e3,
    )
    # ------------------------------------------------------------------- #
    # Physical Parameters ----------------------------------------------- #
    P       = ( 
        k       =   3,              #   Thermal Conductivity [ W/m/K ]
        cp      =   1000,           #   Specific Heat Capacity [ J/kg/K ]
        ρ       =   3200,           #   Density [ kg/m^3 ]
        K0      =   273.15,         #   Kelvin at 0 °C
        Q0      =   0               #   Heat production rate
    )
    P1      = (
        κ       =   P.k/P.ρ/P.cp,   #   Thermal Diffusivity [ m^2/s ] 
    )
    P       =   merge(P,P1)
    # Numerische Konstanten ============================================= #
    NC  =   (
        x       =   100,        # Number of horizontal centroids
        y       =   100,        # Number of vertical centroids
    )
    NV  =   (
        x       =   NC.x + 1,  # Number of horizontal vertices
        y       =   NC.y + 1,  # Number of vertical vertices
    )
    Δ   =   (
        x   =   (abs(M.xmin)+M.xmax)/NC.x,
        y   =   (abs(M.ymin)+M.ymax)/NC.y,
    )
    # ------------------------------------------------------------------- #
    # Erstellung des Gitters ============================================ #
    x   =   (
        c       =   LinRange(M.xmin + Δ.x/2.0, M.xmax - Δ.x/2.0, NC.x),
        ce      =   LinRange(M.xmin - Δ.x/2.0, M.xmax + Δ.x/2.0, NC.x+2),
        v       =   LinRange(M.xmin, M.xmax , NV.x)
    )
    y       = (
        c       =   LinRange(M.ymin + Δ.y/2.0, M.ymax - Δ.y/2.0, NC.y),
        ce      =   LinRange(M.ymin - Δ.x/2.0, M.ymax + Δ.x/2.0, NC.y+2),
        v       =   LinRange(M.ymin, M.ymax, NV.y),
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
    # Tracer Advektions Verfahren ======================================== #
    Ma  =   (
        nmx     =   5,
        nmz     =   5
    )
    # -------------------------------------------------------------------- #
    # Animationssettings ================================================= #
    path        =   string("./exercises/Results/")
    anim        =   Plots.Animation(path, String[] )
    filename    =   string("EnergyEquation",Ini.T,"_",Ini.V,
                            "_",FD.Method.Adv,"_",FD.Method.Diff)
    save_fig    =   0
    # -------------------------------------------------------------------- #
    # Anfangsbedingungen ================================================= #
    # Temperatur --------------------------------------------------------- #
    D       =   (
        Q       =   zeros(NC...),
        T       =   zeros(NC...),
        T0      =   zeros(NC...),
        T_ex    =   zeros(NC.x+2,NC.y+2),
        T_exo   =   zeros(NC.x+2,NC.y+2),
        ρ       =   zeros(NC...),
        cp      =   zeros(NC...),
        vx      =   zeros(NV.x,NV.y+1),
        vy      =   zeros(NV.x+1,NV.y),    
        vxc     =   zeros(NC...),
        vyc     =   zeros(NC...),
        vc      =   zeros(NC...),
        Tmax    =   [0.0],
        Tmin    =   [0.0],
        Tmean   =   [0.0],
    )
    IniTemperature!(Ini.T,M,NC,Δ,D,x,y)
    if FD.Method.Adv==:slf
        D.T_exo    .=  D.T_ex
    end
    # Heat production rate ---
    @. D.Q          = P.Q0

    # Geschwindigkeit ---------------------------------------------------- #
    IniVelocity!(Ini.V,D,NV,Δ,M,x,y)
    # Get the velocity on the centroids ---
    for i = 1:NC.x, j = 1:NC.y
        D.vxc[i,j]  = (D.vx[i,j+1] + D.vx[i+1,j+1])/2
        D.vyc[i,j]  = (D.vy[i+1,j] + D.vy[i+1,j+1])/2
    end
    @. D.vc        = sqrt(D.vxc^2 + D.vyc^2)
    # Visualize initial condition ---------------------------------------- #
    p = heatmap(x.c./1e3 , y.c./1e3, (D.T./D.Tmax)', 
        color=:thermal, colorbar=true, aspect_ratio=:equal, 
        xlabel="x [km]", ylabel="z[km]", 
        title="Temperature", 
        xlims=(M.xmin./1e3, M.xmax./1e3), ylims=(M.ymin./1e3, M.ymax./1e3), 
        clims=(0.5, 1.0))
    quiver!(p,x.c2d[1:Pl.inc:end,1:Pl.inc:end]./1e3,
            y.c2d[1:Pl.inc:end,1:Pl.inc:end]./1e3,
            quiver=(D.vxc[1:Pl.inc:end,1:Pl.inc:end].*Pl.sc,
                D.vyc[1:Pl.inc:end,1:Pl.inc:end].*Pl.sc),        
            color="white")
    if save_fig == 1
        Plots.frame(anim)
    elseif save_fig == 0
        display(p)
    end
    # ------------------------------------------------------------------- #
    # Boundary Conditions ----------------------------------------------- #
    BC     = (type    = (W=:Dirichlet, E=:Dirichlet, 
                        N=:Dirichlet, S=:Dirichlet),
            val     = (W=D.T[1,:],E=D.T[end,:],
                        N=D.T[:,end],S=D.T[:,1]))
    # ------------------------------------------------------------------- #
    # Linear Equations ================================================== #
    if FD.Method.Diff==:implicit || FD.Method.Diff==:CNA
        Num     =   (T=reshape(1:NC.x*NC.y, NC.x, NC.y),)
        ndof    =   maximum(Num.T)        
        if FD.Method.Diff==:CNA
            K1      =   ExtendableSparseMatrix(ndof,ndof)
            K2      =   ExtendableSparseMatrix(ndof,ndof)
        else
            K       =   ExtendableSparseMatrix(ndof,ndof)
        end
        rhs     =   zeros(ndof)
    elseif FD.Method.Diff==:dc
        niter       =   10
        ϵ           =   1e-10
        @. D.ρ      =   P.ρ
        @. D.cp     =   P.cp
        k           =   (x=zeros(NC.x+1,NC.y), y=zeros(NC.x,NC.y+1))
        @. k.x      =   P.k
        @. k.y      =   P.k
        Num         =   (T=reshape(1:NC.x*NC.y, NC.x, NC.y),)
        ndof        =   maximum(Num.T)
        K           =   ExtendableSparseMatrix(ndof,ndof)
        R           =   zeros(NC.x,NC.y)
        ∂T          =   (∂x=zeros(NC.x+1, NC.y), ∂y=zeros(NC.x, NC.y+1))
        q           =   (x=zeros(NC.x+1, NC.y), y=zeros(NC.x, NC.y+1))
    end
    # ------------------------------------------------------------------- #
    # Zeit Konstanten =================================================== #
    T   =   ( 
        tmax    =   [16.336],    # Zeit in Ma
        Δfacc   =   1.0,        # Courant time factor, i.e. dtfac*dt_courant
        Δfacd   =   1.0,        # Diffusion time factor, i.e. dtfac*dt_diff
        Δ       =   [0.0],      # Absolute time step
        Δc      =   [0.0],      # Courant time step
        Δd      =   [0.0],      # Diffusion time stability criterion
    )
    # Courant Kriterium ---
    T.Δc[1]     =   T.Δfacc * minimum((Δ.x,Δ.y)) / 
                        (sqrt(maximum(abs.(D.vx))^2 + maximum(abs.(D.vy))^2))
    T.Δd[1]     =   T.Δfacd * (1.0 / (2.0 * P.κ *(1.0/Δ.x^2 + 1/Δ.y^2)))
    T.Δ[1]      =   minimum([T.Δd[1],T.Δc[1]])
    T.tmax[1]   =   T.tmax[1]*(1e6*(60*60*24*365.25))   # Zeit in [s]
    nt          =   ceil(Int,T.tmax[1]/T.Δ[1])     # Anzahl der Zeitschritte
    # ------------------------------------------------------------------- #
    # Zeitschleife ====================================================== #
    for i=2:nt
        display(string("Time step: ",i))    
        
        # Loesung Advektionsgleichung ---
        if FD.Method.Adv==:upwind
            upwindc2D!(D,NC,T,Δ)
        elseif FD.Method.Adv==:slf
            slfc2D!(D,NC,T,Δ)       
        elseif FD.Method.Adv==:semilag
            semilagc2D!(D,[],[],x,y,T)       
        end

        # Loesung Diffusionsgleichung ---
        if FD.Method.Diff==:explicit
            ForwardEuler2Dc!(D, P.κ, Δ.x, Δ.y, T.Δ[1], P.ρ, P.cp, NC, BC)
        elseif FD.Method.Diff==:implicit
            BackwardEuler2Dc!(D, P.κ, Δ.x, Δ.y, T.Δ[1], P.ρ, P.cp, NC, BC, rhs, K, Num)
        elseif FD.Method.Diff==:CNA
            CNA2Dc!(D, P.κ, Δ.x, Δ.y, T.Δ[1], P.ρ, P.cp, NC, BC, rhs, K1, K2, Num)
        elseif FD.Method.Diff==:ADI
            ADI2Dc!(D, P.κ, Δ.x, Δ.y, T.Δ[1], P.ρ, P.cp, NC, BC)
        elseif FD.Method.Diff==:dc
            D.T0    .=  D.T
            for iter = 1:niter
                # Evaluate residual
                ComputeResiduals2D!(R, D.T, D.T_ex, D.T0, ∂T, q, D.ρ, D.cp, k, BC, Δ, T.Δ[1])
                # @printf("||R|| = %1.4e\n", norm(R)/length(R))
                norm(R)/length(R) < ϵ ? break : nothing
                # Assemble linear system
                K  = AssembleMatrix2D(D.ρ, D.cp, k, BC, Num, NC, Δ, T.Δ[1])
                # Solve for temperature correction: Cholesky factorisation
                Kc = cholesky(K.cscmatrix)
                # Solve for temperature correction: Back substitutions
                δT = -(Kc\R[:])
                # Update temperature
                @. D.T += δT[Num.T]
            end
            D.T0    .= D.T
        end

        display(string("ΔT = ",((D.Tmax[1]-maximum(D.T))/D.Tmax[1])*100))

        # Update Temperatur

        # Plot Solution ---
        if mod(i,10) == 0 || i == nt
            p = heatmap(x.c./1e3 , y.c./1e3, (D.T./D.Tmax)', 
                color=:thermal, colorbar=true, aspect_ratio=:equal, 
                xlabel="x [km]", ylabel="z[km]", 
                title="Temperature", 
                xlims=(M.xmin./1e3, M.xmax./1e3), ylims=(M.ymin./1e3, M.ymax./1e3), 
                clims=(0.5, 1.0))
            quiver!(p,x.c2d[1:Pl.inc:end,1:Pl.inc:end]./1e3,
                    y.c2d[1:Pl.inc:end,1:Pl.inc:end]./1e3,
                    quiver=(D.vxc[1:Pl.inc:end,1:Pl.inc:end].*Pl.sc,
                            D.vyc[1:Pl.inc:end,1:Pl.inc:end].*Pl.sc),        
                color="white")
            if save_fig == 1
                Plots.frame(anim)
            elseif save_fig == 0
                display(p)                        
            end
        end

    end # End Time loop ---

    # Save Animation ---------------------------------------------------- #
    if save_fig == 1
        # Write the frames to a GIF file
        Plots.gif(anim, string( path, filename, ".gif" ), fps = 15)
        foreach(rm, filter(startswith(string(path,"00")), readdir(path,join=true)))
    elseif save_fig == 0
        display(plot(p))
    end

end # Function end ---

EnergyEquation()
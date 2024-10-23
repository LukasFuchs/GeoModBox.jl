using Plots, SpecialFunctions
using GeoModBox.HeatEquation.OneD

function ContinentalGeotherm_1D(T,M,N,Py,t,plotparam)
    # Function to calculate the 1D geotherm for a continental lithosphere.    #
    # Temperature is calculated by solving the 1-D heat equation assuming     #
    # variable thermal parameters and a radiogenic heat source.               #
    # The equation is solved using a proper conserving finite difference      #
    # scheme.                                                                 #
    # The vertical axis is pointing downwards in negative direction.          #
    #                                                                         #
    # Input:                                                                  #
    #   A list of structural parameters M,N,Py,t,plotparam                    #
    #   M - contains:                                                         #
    #       H   -   Model depth [ m ]                                         #
    #       zUC -   Lower upper crust boundary [ m ]                          #
    #       zLC -   Lower lower crust boundary [ m ]                          #
    #   N - contains:                                                         #
    #       nz  -   Number of grid points                                     #
    #   Py - contains:                                                        #
    #       for each layer (upper[UC] and lower[LC] crust and mantle[M])      #
    #       ρ -   Density [ kg/m^3 ]                                        #
    #       cp  -   Specific heat [ J/kg/K ]                                  #
    #       k   -   Thermal conductivity [ W/m/K ]                            #
    #       H   -   Radiogentic heat source per mass [ W/kg ]                 #
    #   t - contains:                                                         #
    #       dtfac   -   Courant criterion                                     #
    #       age     -   Lithospheric age [ Ma ]                               #
    #   plotparam:                                                            #
    #       0   -   no plotting                                               #
    #       1   -   Plot initial and final temperature profile and thermal    #
    #               parameters. The solution of the 1-D heat equation is      #
    #               compared to the steady state solution (poisson equation). #
    #                                                                         #
    # Output:                                                                 #
    #   T - contains:                                                         #
    #       Tpot    -   Potential mantle temperature [ K ]                    #
    #       dTadi   -   Adiabatic mantle temperature gradient [ K/km ]        #
    #       T0      -   Surface temperature [ K ]                             #
    #       T1      -   Bottom temperature [ K ]                              #
    #       T       -   Tempeture profile [ K ]                               #
    #                                                                         #
    # ----------------------------------------------------------------------- #
    #    LF - 17.02.2023 -                                                    #
    # ======================================================================= #
    
    if isempty(T) && isempty(M) && isempty(N) && isempty(Py) && isempty(t)
        ## ================================================================== #
        # Use some default values if no input parameters are defined ======== #
        # Constants --------------------------------------------------------- #
        H           =   200e3               #   Hight of the model [ m ]
        yUC         =   10e3                #   Depth of the upper crust [ m ]
        yLC         =   35e3                #   Depth of the lower crust [ m ]
        
        nc          =   200                 #   Number of grid points
        Δy          =   H/nc                #   Grid resolution
        
        # Depth [m] ---
        yc          =   LinRange(-H + Δy/2.0,0.0 - Δy/2.0,nc)     
        yv          =   LinRange(-H,0,nc+1)
        
        Py  =   (
            # Mantle properties
            ρM      =   3000,               #   Density [ kg/m^3 ]
            cpM     =   1000,               #   Heat capacity [ J/kg/K ]
            kM      =   2.3,                #   Conductivity [ W/m/K ]
            HM      =   0,                  #   Heat generation rate [W/kg]; Q = ρ*H;2.3e-12
        
            # Upper crust properties
            ρUC     =   2700,               # [ kg/m^3 ]
            kUC     =   3.0,                # [ W/m/K ]
            HUC     =   617e-12,            # [ W/kg ]
            cpUC    =   1000,
        
            # Lower crust properties
            ρLC     =   2900,               # [ kg/m^3 ]
            kLC     =   2.0,                # [ W/m/K ]
            HLC     =   43e-12,             # [ W/kg ]
            cpLC    =   1000,
        )                
        # ------------------------------------------------------------------- #
        # Initial Condition ------------------------------------------------- #
        T   =   (
            Tpot    =   1315 + 273.15,      #   Potential temperautre [ K ]
            ΔTadi   =   0.5,                #   Adiabatic temperature gradient [ K/km ]
            Ttop    =   273.15,             #   Surface temperature [ K ]
            T_ex    =   zeros(nc+2,1),    
        )
        T1  =   (
            Tbot    =   T.Tpot + T.ΔTadi*abs(H/1e3),    # Bottom temperature [ K ]
            T       =   T.Tpot .+ abs.(yc./1e3).*T.ΔTadi,   # Initial T-profile [ K ]
        )
        T   =   merge(T,T1)
        
        Tini        =   zeros(nc,1)
        Tini        .=   T.T
        # ------------------------------------------------------------------- #
        # Boundary conditions ----------------------------------------------- #
        BC      =   (
            type    = (N=:Dirichlet, S=:Dirichlet),
            val     = (N=T.Ttop[1],S=T.Tbot[1])
        )
        # S      =   -0.03;          # c     =   -k/q -> 90 mW/m^2
        # N      =   -0.0033;        # c     =   -k/q -> 10 mW/m^2
        # ------------------------------------------------------------------- #
        # =================================================================== #
        # Time stability criterion ------------------------------------------ #
        fac     =   0.8                 #   Courant criterion
        tmax    =   60                  #   Lithosphere age [ Ma ]
        tsca    =   60*60*24*365.25     #   Seconds per year

        age     =   tmax*1e6*tsca        #   Age in seconds    
        # =================================================================== #        
        # Plot Initial condition -------------------------------------------- #
        plotparam   =   1
        p = plot(Tini.-T.Ttop,yc./1e3, 
            label="", 
            xlabel="T [°C]", ylabel="z [m]", 
            title="Initial Temperature",
            xlim=(0,T.Tbot-T.Ttop),ylim=(-H/1e3,0))
        display(p)
        # =================================================================== #
    else
        # =================================================================== #
        # Initial Condition ------------------------------------------------- #
        T1  =   (
            Tbot    =   T.Tpot + T.ΔTadi*abs(H./1e3),    # Bottom temperature [ K ]
            T       =   T.Tpot .+ abs(yc./1e3).*T.ΔTadi,   # Initial T-profile [ K ]
        )
        T   =   merge(T,T1)
        Tini        =   zeros(nc,1)
        Tini        .=   T.T
        # =================================================================== #
        if plotparam==1
            p = plot(T.T.-T.Ttop,yc./1e3, 
            label="", 
            xlabel="T [°C]", ylabel="z [m]", 
            title="Initial Temperature",
            xlim=(0,T.Tbot-T.Ttop),ylim=(-H/1e3,0))
            display(p)
        end
    end
    ## Setup fields ========================================================= #
    Py1     =   (
        k            =   zeros(nc+1,1),
        ρ            =   zeros(nc,1),
        cp           =   zeros(nc,1),
        H            =   zeros(nc,1),
    )
    Py  = merge(Py,Py1)
    
    for j = 1:nc
        if yc[j] >= -yUC
            # Upper Crust ---
            Py.ρ[j]     =   Py.ρUC
            Py.cp[j]    =   Py.cpUC
            Py.H[j]     =   Py.HUC
        elseif yc[j] >= -yLC && yc[j] < -yUC
            # Lower Crust ---
            Py.ρ[j]     =   Py.ρLC
            Py.cp[j]    =   Py.cpLC
            Py.H[j]     =   Py.HLC
        else
            # Mantle ---
            Py.ρ[j]     =   Py.ρM
            Py.cp[j]    =   Py.cpM
            Py.H[j]     =   Py.HM
        end     
    end
    
    for j = 1:nc+1
        if yv[j] >= -yUC
            # Upper Crust ---
            Py.k[j]     =   Py.kLC
        elseif yv[j] >= -yLC && yv[j] < -yUC
            # Lower Crust ---
            Py.k[j]     =   Py.kUC
        else
            # Mantle ---
            Py.k[j]     =   Py.kM
        end 
    end              
    Py2     =   (
        # Thermal diffusivity [ m^2/s ] 
        κ       =  maximum(Py.k)/minimum(Py.ρ)/minimum(Py.cp),     
    )
    Py  =   merge(Py,Py2)
    T2  =   (
        q   =   zeros(nc+1,1),
    )
    T   =   merge(T,T2)  
    # ======================================================================= #
    # Time stability criterion ==0=========================================== #
    Δtexp   =   Δy^2/2/Py.κ             #   Stability criterion for explicit
    Δt      =   fac*Δtexp               #   total time step

    nit     =   ceil(Int64,age/Δt)      #   Number of iterations    

    time    =   zeros(1,nit)            #   Time array
    # ======================================================================= #
    # Calculate 1-D temperature profile ===================================== #
    #count   =   1
    for i = 1:nit
        if i > 1
            time[i]     =   time[i-1] + Δt
        elseif time[i] > age
            Δt          =   age - time[i-1]
            time[i]     =   time[i-1] + Δt
        end
        ForwardEuler1D!(T,Py,Δt,Δy,nc,BC)
        #if i == nit || abs(time[i]/1e6/tsca - count*5.0) < Δt/1e6/tsca        
        #    println(string("i = ",i,", time = ", time[i]/1e6/tsca))
        #    #println(abs(time[i]/1e6/tsca - count*5.0))                
        #    p = plot(T.T.-T.Ttop,yc./1e3, 
        #        label=string("t = ",ceil(time[i]/1e6/tsca),"[Ma]"), 
        #        xlabel="T [°C]", ylabel="z [m]", 
        #        title="Initial Temperature",
        #        xlim=(0,T.Tbot-T.Ttop),ylim=(-H/1e3,0))
        #    display(p)
        #    count = count + 1
        #end
    end
    # ======================================================================= #
    # Calculate heaf flow =================================================== #
    # South ---
    T.T_ex[2:end-1]     =   T.T
    T.T_ex[1]   =   (BC.type.S==:Dirichlet) * (2 * BC.val.S - T.T_ex[2]) + 
                    (BC.type.S==:Neumann) * (T.T_ex[2] - BC.val.S*Δy)
    # North ---
    T.T_ex[end] =   (BC.type.N==:Dirichlet) * (2 * BC.val.N - T.T_ex[nc+1]) +
                    (BC.type.N==:Neumann) * (T.T_ex[nc+1] + BC.val.N*Δy)
    if size(Py.k,1)==1
        for j=1:nc+1
            T.q[j]  =   -Py.k * 
                (T.T_ex[j+1] - T.T_ex[j])/Δy
        end    
    else
        for j=1:nc+1
            T.q[j]  =   -Py.k[j] * 
                (T.T_ex[j+1] - T.T_ex[j])/Δy
        end
    end

    # ======================================================================= #
    ## Plot profile if requested ============================================ #
    if plotparam == 1        
        q = plot(T.T.-T.Ttop,yc./1e3, 
                label=string("t = ",ceil(maximum(time)/1e6/tsca),"[Ma]"), 
                xlabel="T [°C]", ylabel="z [m]",
                title="Initial Temperature",
                xlim=(0,T.Tbot-T.Ttop),ylim=(-H/1e3,0),
                layout=(1,2),subplot=1)        
        q = plot!(T.q.*1e3,yv./1e3, 
                label="", 
                xlabel="q [ mW ]", ylabel="z [m]", 
                title="Heat Flux",
                ylim=(-H/1e3,0),
                subplot=2)        
        display(q)
    end
    # ======================================================================= #
    # keyboard
end
    
T   = []
M   = []
N   = [] 
Py  = []
t   = []
plotparam = 1

ContinentalGeotherm_1D(T,M,N,Py,t,plotparam)
    
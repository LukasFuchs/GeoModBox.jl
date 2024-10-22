using 
function OceanicGeotherm_1D(T,M,N,Py,t,plotparam=0)
# Function to calculate the 1D geotherm for an oceanic lithosphere.       #
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
#   N - contains:                                                         #
#       nz  -   Number of grid points                                     #
#   Py - contains:                                                        #
#       rho -   Density [ kg/m^3 ]                                        #
#       cp  -   Specific heat [ J/kg/K ]                                  #
#       k   -   Thermal conductivity [ W/m/K ]                            #
#       H   -   Radiogenic heat source per mass [ W/kg ]                  #
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
#    LF - 22.10.2024 - juila                                              #
# ======================================================================= #
    
if isempty(M) && isempty(N) && isempty(Py) && isempty(t)
    ## ================================================================== #
    # Use some default values if no input parameters are defined ======== #
    # Constants --------------------------------------------------------- #
    M   =   (
        H       =   -200e3,             # Hight of the model [ m ]
    )
    NC  =   (
        y       =   200,                # Number of grid points
    )
    Δ   =   (
        y       =   M.H/NC.y,            # Grid resolution
    )
    yc          =   LinRange(-M.H + Δ.y/2.0,0.0 - Δ.y/2.0,NC.y)     # Depth [ m ]
    
    Py  =   (
        ρm      =   3000,               # Density [ kg/m^3 ]
        cpm     =   1000,               # Heat capacity [ J/kg/K ]
        km      =   3.0,                # Conductivity [ W/m/K ]
        HM      =   0,                  # Heat generation rate [W/kg]; Q = ρ*H0
    )    
    # Initial Condition ------------------------------------------------- #
    T   =   (
        Tpot    =   1315 + 273.15,      # Potential temperautre [ K ]
        ΔTadi   =   0.5,                # Adiabatic temperature gradient [ K/km ]
        Ttop    =   273.15              # Surface temperature [ K ]
    )
    T1  =   (
        Tbot    =   T.Tpot + T.ΔTadi*abs(M.H/1e3),    # Bottom temperature [ K ]
        T       =   T.Tpot .+ abs.(yc./1e3).*T.ΔTadi,   # Initial T-profile [ K ]
    )
    T   =   merge(T,T1)
      
    Tini        =   T.T

    # Boundary conditions ----------------------------------------------- #
    BC      =   (
        type    = (N=:Dirichlet, S=:Dirichlet),
        val     = (N=:T.Ttop,S=:T.Tbot)
    )
    #T.utbf      =   -0.03;          # c     =   -k/q -> 90 mW/m^2
    #T.ltbf      =   -0.0033;        # c     =   -k/q -> 10 mW/m^2
    # ------------------------------------------------------------------- #    
    # Time stability criterion ------------------------------------------ #
    t       =   (
        fac     =   0.8,                #   Courant criterion
        max     =   60,                 #   Lithosphere age [ Ma ]
        sca     =   60*60*24*365.25,    #   Seconds per year
    )
    t1      =   (
        age     =   t.max*1e6*t.sca,    #   Age in seconds
    )
    t   =   merge(t,t1)
    
    # =================================================================== #    
    
    ## Plot Initial condition ------------------------------------------- #
    plotparam   =   1
    #fig = figure;
    #clf
    #subplot(1,2,1)
    #plot(T.T,M.z./1e3,'LineWidth',2)
    #xlabel('$$T\ [K]$$','Interpreter','latex')
    #ylabel('$$Depth\ [km]$$','Interpreter','latex')
    #title('$$T-profile$$','Interpreter','latex')
    #set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',15,...
    #    'TickLabelInterpreter','latex')
    # =================================================================== #
else
    # =================================================================== #
    # Initial Condition ------------------------------------------------- #    
    T1  =   (
        Tbot    =   T.Tpot + T.ΔTadi*abs(M.H./1e3),    # Bottom temperature [ K ]
        T       =   T.Tpot .+ abs(yc./1e3).*T.ΔTadi,   # Initial T-profile [ K ]
    )
    T   =   merge(T,T1)
      
    Tini        =   T.T
    # =================================================================== #
    #if plotparam
    #    fig = figure;
    #    clf
    #    subplot(1,2,1)
    #    plot(T.T,M.z./1e3,'LineWidth',2)
    #    xlabel('$$T\ [K]$$','Interpreter','latex')
    #    ylabel('$$Depth\ [km]$$','Interpreter','latex')
    #    title('$$T-profile$$','Interpreter','latex')
    #    set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',15,...
    #        'TickLabelInterpreter','latex')
    #end
end
# Setup Fields ========================================================= #
Py1     =   (
    ρ       =   Py.ρm.*ones(NC.y,1),
    cp      =   Py.cpm.*ones(NC.y,1),
    k       =   Py.km.*ones(NC.y,1),
    H       =   Py.HM.*one(NC.y,1)
)
Py  =   merge(Py,Py1)
Py2     =   (
    κ       =  maximum(Py.k)/minimum(Py.rho)/minimum(Py.cp)     # Therma  l diffusivity [ m^2/s ] 
)
Py  =   merge(Py,Py2)
T2  =   (
    q   =   zeros(NC.y,1)
)
T   =   merge(T,T2)
# ======================================================================= #
## Time stability criterion ============================================= #
t2  =   (
    Δexp    =   N.dz^2/2/Py.kappa;  #   Stability criterion for explicit
)
t   =   merge(t,t2)
t3  =   (
    Δ   =   t.fac*t.Δexp,
    nit =   ceil(t.age./t.dt),      #   Number of iterations
)
t   =   merge(t,t3)
t4  =   (
    time    =   zeros(1,t.nit),     #   Time array
)
t   =   merge(t,t4)
# ======================================================================= #
# Calculate 1-D temperature profile ===================================== #
for i = 1:t.nit
    if i > 1
        t.time[i]   =   t.time[i-1] + t.dt
    end
    println(t.time)
    #     [T]     =   SolveDiff1Dimplicit_vary(N,T,Py,t);
    # [T]     =   SolveDiff1Dexplicit_vary(N,T,Py,t);
end
# ======================================================================= #
# Calculate heaf flow =================================================== #
#if size(Py.k,1)==1
#    for j=1:NC.y
#        T.q[j] = -Py.k * 
#            (T.T[j+1] - T.T[j])/N.dz
#    end
#    T.q[1]          =   - Py.k*(T.T[2]-T.T[1])/Δ.y
#    T.q[NC.y]       =   - Py.k*(T.T[NC.y]-T.T[N.nz-1])/Δ.y
#else
    #for j=1:NC.y
    #    T.q(j) = -(Py.k(j+1) + Py.k(j))/2 * 
    #        (T.T(j+1) - T.T(j))/N.dz;
    #end
    #T.q(1)          =   - Py.k(1)*(T.T(2)-T.T(1))/N.dz;
    #T.q(N.nz-1)     =   - Py.k(N.nz)*(T.T(N.nz)-T.T(N.nz-1))/N.dz;
#end
# ======================================================================= #
    
### Plot profile if requested ============================================ #
#if plotparam
#    if strcmp(T.lbound,'const') && strcmp(T.ubound,'const')
#        T.Tana      =   zeros(N.nz,1);
#        T.Tana      =   T.Tini + ...
#            (T.T0 - T.Tpot).*erfc(-M.z./(2*sqrt(t.age*Py.kappa)));
#        T.Tana(1)   =   T.T0;
#    end
#    figure(fig)
#    subplot(1,2,1)
#    hold on
#    plot(T.T,M.z./1e3,'r-','LineWidth',2)
#    if strcmp(T.lbound,'const') && strcmp(T.ubound,'const')
#        plot(T.Tana,M.z/1e3,'--','LineWidth',2)
#        legend('$$Initial$$',['$$T_{',num2str(t.age/1e6/t.tfac),'Ma}$$'],...
#            '$$T_{HSCM}$$','Location','SouthWest','Interpreter','latex')
#    else
#        legend('$$Initial$$',['$$T_{',num2str(t.age/1e6/t.tfac),'Ma}$$'],...
#            'Location','SouthWest','Interpreter','latex')
#    end
#    xlabel('$$T\ [ K ]$$','Interpreter','latex')
#    ylabel('$$Depth\ [ km ]$$','Interpreter','latex')
#    subplot(1,2,2)
#    plot(T.q.*1e3,M.zc./1e3,'LineWidth',2)
#    xlabel('$$q\ [ mW ]$$','Interpreter','latex')
#    ylabel('$$Depth\ [ km ]$$','Interpreter','latex')
#    title('$$Heat\ flux$$','Interpreter','latex')
#    set(gca,'FontWeight','Bold','LineWidth',2,'FontSize',15,...
#            'TickLabelInterpreter','latex')
#end
# ======================================================================= #
# keyboard
end
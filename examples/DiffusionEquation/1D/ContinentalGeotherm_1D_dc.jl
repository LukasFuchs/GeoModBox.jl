using Plots, ExtendableSparse, LinearAlgebra
using GeoModBox.HeatEquation.OneD
using TimerOutputs
#CHECK HEAT FLUX CALUCLATION!!!

function ContinentalGeotherm_1D_dc()
to = TimerOutput()
# ------------------------------------------------------------------- #
#    LF - 18.12.2025 - julia                                          #
# ------------------------------------------------------------------- #
# Constants --------------------------------------------------------- #
H           =   200e3               #   Hight of the model [ m ]
xUC         =   10e3                #   Depth of the upper crust [ m ]
xLC         =   35e3                #   Depth of the lower crust [ m ]
nc          =   200                 #   Number of grid points
Δx          =   H/nc                #   Grid resolution
# Depth [m] ---
xc          =   LinRange(-H + Δx/2.0,0.0 - Δx/2.0,nc)     
xv          =   LinRange(-H,0,nc+1)
Py  =   (
    # Mantle properties ---
    ρM      =   3000,               #   Density [ kg/m^3 ]
    cpM     =   1000,               #   Heat capacity [ J/kg/K ]
    kM      =   2.3,                #   Conductivity [ W/m/K ]
    HM      =   2.3e-12,            #   Heat generation rate [W/kg]; Q = ρ*H;2.3e-12
        
    # Upper crust properties ---
    ρUC     =   2700,               # [ kg/m^3 ]
    kUC     =   3.0,                # [ W/m/K ]
    HUC     =   617e-12,            # [ W/kg ]
    cpUC    =   1000,
        
    # Lower crust properties ---
    ρLC     =   2900,               # [ kg/m^3 ]
    kLC     =   2.0,                # [ W/m/K ]
    HLC     =   43e-12,             # [ W/kg ]
    cpLC    =   1000,
)                
# Iterations ---
niter       =   50  
ϵ           =   1.0e-15 
R           =   zeros(nc)   
# Discretization method --- 
C           =   1.0
# Assemble Coefficient Matrix --------------------------------------- #
ndof        =   nc
K           =   ExtendableSparseMatrix(ndof,ndof)    
# ------------------------------------------------------------------- #     
# Initial Condition ------------------------------------------------- #
T   =   (
    Tpot    =   1315 + 273.15,      #   Potential temperautre [ K ]
    ΔTadi   =   0.5,                #   Adiabatic temperature gradient [ K/km ]
    Ttop    =   273.15,             #   Surface temperature [ K ]
    T_ex    =   zeros(nc+2),    
    T0      =   zeros(nc),
    T_ex0   =   zeros(nc+2),
)
T1  =   (
    Tbot    =   T.Tpot + T.ΔTadi*abs(H/1e3),    # Bottom temperature [ K ]
    T       =   T.Tpot .+ abs.(xc./1e3).*T.ΔTadi,   # Initial T-profile [ K ]
)
T               =   merge(T,T1)
Tini            =   zeros(nc)
Tini            .=   T.T
T.T_ex[2:end-1] .=   T.T
T.T0            .=  T.T
T.T_ex0         .=  T.T_ex
# ------------------------------------------------------------------- #
# Boundary conditions ----------------------------------------------- #
BC      =   (
    type    = (W=:Dirichlet, E=:Dirichlet),
    val     = (W=T.Tbot,E=T.Ttop)
)
# S      =   -0.03;          # c     =   -k/q -> 90 mW/m^2
# N      =   -0.0033;        # c     =   -k/q -> 10 mW/m^2
# ------------------------------------------------------------------- #
# Time stability criterion ------------------------------------------ #
fac     =   0.9                 #   Courant criterion
tmax    =   1000                 #   Lithosphere age [ Ma ]
tsca    =   60*60*24*365.25     #   Seconds per year
age     =   tmax*1e6*tsca        #   Age in seconds    
# ------------------------------------------------------------------- #
# Plot Initial condition -------------------------------------------- #
p = plot(T.T,xc./1e3, 
    label="", 
    xlabel="T [ K ]", ylabel="y [ km ]", 
    title="Initial Temperature",
    xlim=(T.Ttop,T.Tbot),ylim=(-H/1e3,0))
display(p)
# ------------------------------------------------------------------- #
# Setup fields ------------------------------------------------------ #
Py1     =   (
    k            =   zeros(nc+1),
    ρ            =   zeros(nc),
    cp           =   zeros(nc),
    H            =   zeros(nc),
    )
Py  = merge(Py,Py1)
# Upper Crust ---
Py.k[xv.>=-xUC]     .=   Py.kUC
Py.ρ[xc.>=-xUC]     .=   Py.ρUC
Py.cp[xc.>=-xUC]    .=   Py.cpUC
Py.H[xc.>=-xUC]     .=   Py.HUC    
# Lower Crust ---
Py.k[xv.>=-xLC .&& xv.<-xUC]     .=   Py.kLC
Py.ρ[xc.>=-xLC .&& xc.<-xUC]     .=   Py.ρLC
Py.cp[xc.>=-xLC .&& xc.<-xUC]    .=   Py.cpLC
Py.H[xc.>=-xLC .&& xc.<-xUC]     .=   Py.HLC   
# Mantle ---
Py.k[xv.<-xLC]  .=   Py.kM
Py.ρ[xc.<-xLC]  .=   Py.ρM
Py.cp[xc.<-xLC] .=   Py.cpM
Py.H[xc.<-xLC]  .=   Py.HM
Py2     =   (
    # Thermal diffusivity [ m^2/s ] 
    κ       =   maximum(Py.k)/minimum(Py.ρ)/minimum(Py.cp),  
    Q       =   zeros(nc),
    )
Py  =   merge(Py,Py2)  
@. Py.Q    =   Py.H * Py.ρ
q   =   (
    x   =   zeros(nc+1),
    x0  =   zeros(nc+1),
)
∂T  =   (
    ∂x  =   zeros(nc+1),
    ∂x0 =   zeros(nc+1),
)
# ------------------------------------------------------------------- #
# Time stability criterion ------------------------------------------ #
Δtexp   =   Δx^2/2/Py.κ             #   Stability criterion for explicit
Δt      =   fac*Δtexp               #   total time step
nit     =   ceil(Int64,age/Δt)      #   Number of iterations    
time    =   zeros(1,nit)            #   Time array
# ------------------------------------------------------------------- #
# Time Loop --------------------------------------------------------- #
count   =   1
@timeit to "TimeLoop" begin
for i = 1:nit
    if i > 1
        time[i]     =   time[i-1] + Δt
    elseif time[i] > age
        Δt          =   age - time[i-1]
        time[i]     =   time[i-1] + Δt
    end
    for iter = 1:niter
        ComputeResiduals1D!(R, T.T, T.T_ex, T.T0, T.T_ex0, Py.Q, ∂T,
                 q, Py.ρ, Py.cp, Py.k, BC, Δx, Δt;C)
        # @printf("||R|| = %1.4e\n", norm(R)/length(R))            
        norm(R)/length(R) < ϵ ? break : nothing
        # Assemble linear system
        AssembleMatrix1D( Py.ρ, Py.cp, Py.k, Δx, Δt, nc, BC, K;C)
        # Solve for temperature correction: Cholesky factorisation
        Kc = cholesky(K.cscmatrix)
        # Solve for temperature correction: Back substitutions
        δT = -(Kc\R[:])          
        # Update temperature            
        T.T .= T.T .+ δT  
    end
    T.T0    .= T.T
    if i == nit || abs(time[i]/1e6/tsca - count*100.0) < Δt/1e6/tsca        
        println(string("i = ",i,", time = ", time[i]/1e6/tsca))                    
        p = plot!(T.T.-T.Ttop,xc./1e3, 
            label=string("t = ",ceil(time[i]/1e6/tsca),"[Ma]"), 
            xlabel="T [°C]", ylabel="z [m]", 
            title="Continental Geotherm",
            xlim=(0,T.Tbot-T.Ttop),ylim=(-H/1e3,0))
        display(p)
        count = count + 1
    end
end
end
# ------------------------------------------------------------------- #
# Plot profile if requested ----------------------------------------- #
p1 = plot(T.T.-T.Ttop,xc./1e3, 
        label=string("t = ",ceil(maximum(time)/1e6/tsca),"[Ma]"), 
        xlabel="T [°C]", ylabel="z [m]",
        title="Continental Geotherm",
        xlim=(0,T.Tbot-T.Ttop),ylim=(-H/1e3,0),
        layout=(1,3),subplot=1)    
p1 = plot!(q.x.*1e3,xv./1e3, 
        label="", 
        xlabel="q [ mW ]", ylabel="z [m]", 
        title="Heat Flux",
        ylim=(-H/1e3,0),
        subplot=2) 
p1 = plot!(Py.k,xv./1e3,label="k [W/m/K]",
        xlabel="k,ρ,cp,Q", ylabel="z [km]",
        title="Thermal parameters",
        subplot=3)       
p1 = plot!(Py.cp./1e3,xc./1e3,label="cp [kJ/kg/K]",
        subplot=3)
p1 = plot!(Py.ρ,xc./1e3,label="ρ [kg/m³]",                
        subplot=3)
p1 = plot!(Py.H.*(Py.ρ./1e3),xc./1e3,label="Q [mW/m³]",
        subplot=3)
display(p1)
# Save figures ---
savefig(p1,"./examples/DiffusionEquation/1D/Results/ContinentalGeotherm_1D_dc.png")
savefig(p,"./examples/DiffusionEquation/1D/Results/ContinentalGeotherm_1D_evolve_dc.png")
# ------------------------------------------------------------------- #
display(to)
end
    
ContinentalGeotherm_1D_dc()
    
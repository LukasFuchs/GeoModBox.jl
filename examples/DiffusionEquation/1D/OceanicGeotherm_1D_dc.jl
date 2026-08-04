using Plots, SpecialFunctions
using GeoModBox.HeatEquation.OneD
using TimerOutputs, ExtendableSparse, Printf, LinearAlgebra

function OceanicGeotherm_1D_dc()
to  =   TimerOutput()
# ------------------------------------------------------------------- #
#    LF - 18.12.2025 - juila                                          #
# ------------------------------------------------------------------- #
# Constants --------------------------------------------------------- #
H           =   200e3               #   Hight of the model [ m ]
nc          =   200                 #   Number of central grid points    
Δx          =   H/nc                #   Grid resolution
# Depth [ m ] ---
xc          =   LinRange(-H + Δx/2.0,0.0 - Δx/2.0,nc)     
xv          =   LinRange(-H,0,nc+1)
# Physical Parameters ---
Py  =   (
    ρm      =   3000,               #   Density [ kg/m^3 ]
    cpm     =   1000,               #   Heat capacity [ J/kg/K ]
    km      =   3.0,                #   Conductivity [ W/m/K ]
    HM      =   0,                  #   Heat generation rate [W/kg]; Q = ρ*H0
)
# Iterations ---
niter       =   50  
ϵ           =   1.0e-15 
R           =   zeros(nc)     
# Discretization Method ---
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
Tini            .=  T.T
T.T_ex[2:end-1] .=  T.T
T.T0            .=  T.T
T.T_ex0         .=  T.T_ex
# ------------------------------------------------------------------- #
# Boundary conditions ----------------------------------------------- #
BC      =   (
    type    = (W=:Dirichlet, E=:Dirichlet),
    val     = (W=T.Tbot[1],E=T.Ttop[1])
)
# S      =   -0.03;          # c     =   -k/q -> 90 mW/m^2
# N      =   -0.0033;        # c     =   -k/q -> 10 mW/m^2
# ------------------------------------------------------------------- #
# Time stability criterion ------------------------------------------ #
fac     =   0.8                 #   Courant criterion
tmax    =   60                  #   Lithosphere age [ Ma ]
tsca    =   60*60*24*365.25     #   Seconds per year
age     =   tmax*1e6*tsca        #   Age in seconds    
# ------------------------------------------------------------------- #
# Plot Initial condition -------------------------------------------- #
p1 = plot(Tini,xc./1e3, 
    label="", 
    xlabel="T [ K ]", ylabel="y [ km ]", 
    title="Initial Temperature",
    xlim=(T.Ttop,T.Tbot),ylim=(-H/1e3,0))
display(p1)
# ------------------------------------------------------------------- #
# Setup Fields ------------------------------------------------------ #
Py1     =   (
    ρ       =   Py.ρm.*ones(nc),
    cp      =   Py.cpm.*ones(nc),
    k       =   Py.km.*ones(nc+1),
    H       =   Py.HM.*ones(nc)
)
Py  =   merge(Py,Py1)
Py2     =   (
    # Thermal diffusivity [ m^2/s ] 
    κ       =  maximum(Py.k)/minimum(Py.ρ)/minimum(Py.cp),     
)
Py  =   merge(Py,Py2)
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
# Calculate 1-D temperature profile --------------------------------- #
@timeit to "TimeLoop" begin
count   =   1
for i = 1:nit
    if i > 1
        time[i]     =   time[i-1] + Δt
    elseif time[i] > age
        Δt          =   age - time[i-1]
        time[i]     =   time[i-1] + Δt
    end
    for iter = 1:niter
        ComputeResiduals1D!(R, T.T, T.T_ex, T.T0, T.T_ex0, Py.H, ∂T,
                 q, Py.ρ, Py.cp, Py.k, BC, Δx, Δt;C)
        @printf("||R|| = %1.4e\n", norm(R)/length(R))            
        norm(R)/length(R) < ϵ ? break : nothing
        # Assemble linear system
        AssembleMatrix1D( Py.ρ, Py.cp, Py.k, Δx, Δt, nc, BC, K;C)
        # Solve for temperature correction: Cholesky factorisation
        Kc = cholesky(K.cscmatrix)
        # Solve for temperature correction: Back substitutions
        δT = -(Kc\R[:])          
        # Update temperature            
        T.T .= T.T .+ δT  
        # ForwardEuler1D!(T,Py,Δt,Δy,nc,BC)
    end
    T.T0    .= T.T
    if i == nit || abs(time[i]/1e6/tsca - count*5.0) < Δt/1e6/tsca        
        println(string("i = ",i,", time = ", time[i]/1e6/tsca))        
        p1 = plot!(T.T,xc./1e3, 
            label=string("t = ",ceil(time[i]/1e6/tsca),"[Ma]"), 
            xlabel="T [ K ]", ylabel="y [ km ]", 
            title="Oceanic Geotherm")#,
            # xlim=(T.Ttop,T.Tbot),ylim=(-H/1e3,0))
        display(p1)
        count = count + 1
    end
end
end
# ------------------------------------------------------------------- #
# Plot -------------------------------------------------------------- #
if BC.type.W==:Dirichlet && BC.type.E==:Dirichlet
    Tana    =   zeros(nc,1)
    Tana    .=   Tini .+ (T.Ttop - T.Tpot).*erfc.(-xc./(2*sqrt(age*Py.κ)))
    Tana[end] =   T.Ttop
end    
p = plot(T.T,xc./1e3, 
        label=string("t = ",ceil(maximum(time)/1e6/tsca),"[Ma]"), 
        xlabel="T [ K ]", ylabel="y [ km ]",
        title="Oceanic Geotherm",
        xlim=(T.Ttop,T.Tbot),ylim=(-H/1e3,0),
        layout=(1,2),subplot=1)        
if BC.type.W==:Dirichlet && BC.type.E==:Dirichlet
    plot!(p,Tana,xc./1e3, 
            label="T_HSCM",linestyle=:dash,
            layout=(1,2),subplot=1)        
end        
p = plot!(q.x.*1e3,xv./1e3, 
        label="", 
        xlabel="q [ mW ]", ylabel="y [m]", 
        title="Heat Flux",
        ylim=(-H/1e3,0),
        subplot=2)        
display(p)
savefig(p,"./examples/DiffusionEquation/1D/Results/OceanicGeotherm_1D_dc.png")
savefig(p1,"./examples/DiffusionEquation/1D/Results/OceanicGeotherm_1D_dc_evolve.png")
# ======================================================================= #
display(to)
end

OceanicGeotherm_1D_dc()
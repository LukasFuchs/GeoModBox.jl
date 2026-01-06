# ======================================================================= #
# Lösen der 1-D Diffusionsgleichung                                       #
# ----------------------------------------------------------------------- #
# LF - 18.12.2025 - Vers. 1.1. - julia                                    #
# ======================================================================= #
using Plots, Printf, LinearAlgebra, ExtendableSparse
using GeoModBox.HeatEquation.OneD
using TimerOutputs

function Heat_1D_dc()
to  =   TimerOutput() 
@timeit to "Ini" begin
# Physical Parameters --------------------------------------------------- #
L           =   100.0               # Length [ m ]
Trock       =   300.0               # Background temperature [ C ]
Tmagma      =   1200.0              # Dike temperature [ C ]
W           =   5.0                 # Dike width [m]
κ           =   1.0e-6              # Diffusivity [ m²/s ]
# ----------------------------------------------------------------------- #
# Numerical Parameters -------------------------------------------------- #
nc          =   100                 # Number of cenroids
Δx          =   L/nc                # Grid spacing
xc          =   Δx/2:Δx:(L-Δx/2)    # Coordinates
# Iterations ---
niter       =   50  
ϵ           =   1.0e-15       
# ----------------------------------------------------------------------- #
# Time Parameters ------------------------------------------------------- #
day         =   3600.0*24.0         # Seconds per day
tmax        =   2.0*365.25*day      # Maximum time [ s ]
time        =   0.0                 # Initilalize time
# Explicit time stability criterion ---
Δtexp       =   Δx^2/κ/2.0
# Absolut time step ---
fac         =   0.9                 # Factorisation
Δt          =   fac*Δtexp           # Absolut time step
# Number of time steps ---
nt          =   ceil(Int,tmax/Δt)
# ----------------------------------------------------------------------- #
# Setting up field memroy ----------------------------------------------- #
∂2T =   (
    ∂x2     = zeros(nc), 
    ∂x20    = zeros(nc), 
)
T   =   (
    ini     = zeros(nc), 
    ana     = zeros(nc),)
expl     = (
    T       = zeros(nc), 
    T0      = zeros(nc), 
    T_ex    = zeros(nc+2),
    T_ex0   = zeros(nc+2),
    ε       = zeros(nc),
    R       = zeros(nc),
)
imp     = (
    T       = zeros(nc), 
    T0      = zeros(nc), 
    T_ex    = zeros(nc+2),
    T_ex0   = zeros(nc+2),
    ε       = zeros(nc),
    R       = zeros(nc),
)
CNA     = (
    T       = zeros(nc), 
    T0      = zeros(nc), 
    T_ex    = zeros(nc+2),
    T_ex0   = zeros(nc+2),
    ε       = zeros(nc),
    R       = zeros(nc)
)
Q   =   zeros(nc)
ρ   =   ones(nc).*3300
cp  =   ones(nc).*1200
# ----------------------------------------------------------------------- #
# Initial condition ----------------------------------------------------- #
# Gaussian temperature distribution ---------
σ           =   5
xp          =   L/2.0
@. T.ini    =   Trock + (Tmagma-Trock)*exp(-((xc-xp)/σ)^2)
expl.T0     .=  T.ini
imp.T0      .=  T.ini
CNA.T0      .=  T.ini
# Analytical solution ---
@. T.ana    =   Trock + (Tmagma-Trock)/(sqrt(1+4*time*κ/σ^2))*
                        exp(-(xc-xp)^2/(σ^2 + 4*time*κ))
# ----------------------------------------------------------------------- #
end
# Boundary conditions --------------------------------------------------- #
BC          =   (
                    type = (W=:Dirichlet, E=:Dirichlet),
                    #type = (W=:Neumann, E=:Neumann),
                    val = (W=:300.0,E=:300.0)
)
# ----------------------------------------------------------------------- #
# Assemble Coefficient Matrix ------------------------------------------- #
ndof        =   length(T.ini)
K           =   ExtendableSparseMatrix(ndof,ndof)    
# ----------------------------------------------------------------------- #
# Animationssettings ---------------------------------------------------- #
path        =   string("./examples/DiffusionEquation/1D/Results/")
anim        =   Plots.Animation(path, String[] )
filename    =   string("1D_comparison_general_solver")
save_fig    =   1
# ----------------------------------------------------------------------- #
# Plot initial condition ------------------------------------------------ #
p = plot(xc, expl.T0, label="explicit", 
        xlabel="x [m]", ylabel="T [°C]", 
        title="Temperature after $(round(time / day, digits=1)) days
        Δt = $(round(Δt / Δtexp, digits=2))*Δt_{crit}",
        xlim=(0,L),ylim=(0, Tmagma),layout=(1,2))
plot!(p,xc, imp.T0,label="implicit",subplot=1)
plot!(p,xc, CNA.T,label="cna",subplot=1)
plot!(p,xc, T.ana, linestyle=:dash, label="analytical",subplot=1)
plot!(p,xc, expl.ε, xlabel="x [m]", ylabel="ε",
        title="Error",
        label="ε_exp",xlim=(0,L),ylim=(0,2.0),
        subplot=2)        
plot!(p,xc, imp.ε, label="ε_imp",subplot=2)      
plot!(p,xc, CNA.ε, label="ε_cna",subplot=2)  
if save_fig == 1
    Plots.frame(anim)
else
    display(p)
end
# ----------------------------------------------------------------------- #
# Time loop ------------------------------------------------------------- #
@timeit to "TimeLoop" begin
for n=1:nt
    println("Zeitschritt: ",n,", Time: $(round(time/day, digits=1)) [d]")
    @timeit to "solution" begin
    for iter = 1:niter
        # Residual iteration
        ComputeResiduals1Dc!( expl.R, expl.T, expl.T_ex, expl.T0, expl.T_ex0, 
                ∂2T, Q, ρ, cp, κ, BC, Δx, Δt; C=1.0 )
        @printf("||R|| = %1.4e\n", norm(expl.R)/length(expl.R))            
        norm(expl.R)/length(expl.R) < ϵ ? break : nothing
        # Assemble linear system
        AssembleMatrix1Dc!( κ, Δx, Δt, nc, BC, K; C=1.0 )
        # Solve for temperature correction: Cholesky factorisation
        Kc = cholesky(K.cscmatrix)
        # Solve for temperature correction: Back substitutions
        δT = -(Kc\expl.R[:])
        # Update temperature            
        expl.T .= expl.T .+ δT            
    end
    # Update temperature ------------------------------------------------ #
    expl.T0    .=  expl.T
    for iter = 1:niter
        # Residual iteration
        ComputeResiduals1Dc!( imp.R, imp.T, imp.T_ex, imp.T0, imp.T_ex0, 
                ∂2T, Q, ρ, cp, κ, BC, Δx, Δt; C=0.0 )
        @printf("||R|| = %1.4e\n", norm(imp.R)/length(imp.R))            
        norm(imp.R)/length(imp.R) < ϵ ? break : nothing
        # Assemble linear system
        AssembleMatrix1Dc!( κ, Δx, Δt, nc, BC, K; C=0.0 )
        # Solve for temperature correction: Cholesky factorisation
        Kc = cholesky(K.cscmatrix)
        # Solve for temperature correction: Back substitutions
        δT = -(Kc\imp.R[:])          
        # Update temperature            
        imp.T .= imp.T .+ δT            
    end
    # Update temperature ------------------------------------------------ #
    imp.T0    .=  imp.T
    for iter = 1:niter
        # Residual iteration
        ComputeResiduals1Dc!( CNA.R, CNA.T, CNA.T_ex, CNA.T0, CNA.T_ex0, 
                ∂2T, Q, ρ, cp, κ, BC, Δx, Δt; C=0.5 )
        @printf("||R|| = %1.4e\n", norm(CNA.R)/length(CNA.R))            
        norm(CNA.R)/length(CNA.R) < ϵ ? break : nothing
        # Assemble linear system
        AssembleMatrix1Dc!( κ, Δx, Δt, nc, BC, K; C=0.5 )
        # Solve for temperature correction: Cholesky factorisation
        Kc = cholesky(K.cscmatrix)
        # Solve for temperature correction: Back substitutions
        δT = -(Kc\CNA.R[:])          
        # Update temperature            
        CNA.T .= CNA.T .+ δT            
    end
    # Update temperature ------------------------------------------------ #
    CNA.T0    .=  CNA.T
    end
    # Update time ------------------------------------------------------- #
    time    =   time + Δt
    # Analytical Solution ----------------------------------------------- #
    @. T.ana    =   Trock + (Tmagma-Trock)/(sqrt(1+4*time*κ/σ^2))*
                        exp(-(xc-xp)^2/(σ^2 + 4*time*κ))
    # Error ------------------------------------------------------------- #
    @. expl.ε    =   abs((T.ana-expl.T0)/T.ana)*100
    @. imp.ε    =   abs((T.ana-imp.T0)/T.ana)*100
    @. CNA.ε    =   abs((T.ana-CNA.T0)/T.ana)*100
    # Plot solution ----------------------------------------------------- #
    if n == 1 || n % 5 == 0 || n == nt
        # Subplot 1 ---
        p = plot(xc, expl.T, label="numerical",
                xlim=(0,L),ylim=(0,1300),
                xlabel="x [m]",ylabel="T [°C]",
                title="Temperature after $(round(time / day, digits=1)) days
                Δt = $(round(Δt / Δtexp, digits=2))*Δt_{crit}",
                layout=(1,2))
        plot!(p, xc, imp.T,linestyle=:dash, label="implicit",subplot=1)
        plot!(p, xc, CNA.T,linestyle=:dash, label="cna",subplot=1)
        plot!(p, xc, T.ana, linestyle=:dash, label="analytical",subplot=1)    
        # Subplot 2 ---
        plot!(p,xc, expl.ε, label="ε_exp",
            xlim=(0,L),ylim=(0,2.0),
            xlabel="x [m]",ylabel="ε [%]",
            title="Error",
            subplot=2)
        plot!(p, xc, imp.ε, label="ε_imp",subplot=2)
        plot!(p, xc, CNA.ε, label="ε_cna",subplot=2)                
        # Display the plots ---    
        if save_fig == 1
            Plots.frame(anim)
        else
            display(p)
        end
    end
end
end
# Speicher Animation ---------------------------------------------------- #
if save_fig == 1
    # Write the frames to a GIF file
    Plots.gif(anim, string( path, filename, ".gif" ), fps = 15)
else
    display(plot(p))
end
foreach(rm, filter(startswith(string(path,"00")), readdir(path,join=true)))
# ----------------------------------------------------------------------- #
display(to)
end
# Call function --------------------------------------------------------- #
Heat_1D_dc()
# ----------------------------------------------------------------------- #
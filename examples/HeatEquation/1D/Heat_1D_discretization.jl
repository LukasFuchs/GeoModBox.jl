# ======================================================================= #
# Lösen der 1-D Diffusionsgleichung                                       #
# ----------------------------------------------------------------------- #
# LF - 18.09.2024 - Vers. 1.0. - julia                                    #
# ======================================================================= #
using Plots, Printf, LinearAlgebra, ExtendableSparse
using GeoModBox.HeatEquation.OneD

function Heat_1D_discretization()
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
# Iterations
niter       =   10  
ϵ           =   1.0e-10       
R           =   zeros(nc)           # Residual
# Time Parameters ------------------------------------------------------- #
day         =   3600.0*24.0         # Seconds per day
tmax        =   2.0*365.25*day      # Maximum time [ s ]
time        =   0.0                 # Initilalize time

Δtexp       =   Δx^2/κ/2.0          # Explicit time stability criterion

fac         =   0.9                 # Factorisation
Δt          =   fac*Δtexp           # Absolut time step

nt          =   ceil(Int,tmax/Δt)   # Number of time steps
# ----------------------------------------------------------------------- #
# Initial condition ----------------------------------------------------- #
T           =   zeros(nc)
T_ex        =   zeros(nc+2)
∂T2∂x2      =   zeros(nc)
# Gaussian temperature distribution ---------
σ           =   5
xp          =   L/2.0
@. T        =   Trock + (Tmagma-Trock)*exp(-((xc-xp)/σ)^2)

T0exp       =   T
T1exp       =   T0exp
T0imp       =   T
T1imp       =   T0imp
T0dc        =   T
T1dc        =   T0dc
T0cnv       =   T
T1cnv       =   T0cnv
Tana        =   zeros(nc)
εexp        =   zeros(nc)
εimp        =   zeros(nc)
εdc         =   zeros(nc)
εcnv        =   zeros(nc)

# Analytical solution ---
@. Tana     =   Trock + (Tmagma-Trock)/(sqrt(1+4*time*κ/σ^2))*
                        exp(-(xc-xp)^2/(σ^2 + 4*time*κ))
# Boundary conditions --------------------------------------------------- #
BC          =   (
                    type = (W=:Dirichlet, E=:Dirichlet),
                    #type = (W=:Neumann, E=:Neumann),
                    val = (W=:300.0,E=:300.0)
)
# ----------------------------------------------------------------------- #
## Assemble Coefficient Matrix ------------------------------------------- #
ndof        =   length(T)
K           =   ExtendableSparseMatrix(ndof,ndof)    
K1          =   ExtendableSparseMatrix(ndof,ndof)    
K2          =   ExtendableSparseMatrix(ndof,ndof)    
## ----------------------------------------------------------------------- #
# Animationssettings ---------------------------------------------------- #
path        =   string("./examples/HeatEquation/1D/Figures/")
anim        =   Plots.Animation(path, String[] )
filename    =   string("1D_comparison")
save_fig    =   1
# ----------------------------------------------------------------------- #
# Plot initial condition ------------------------------------------------ #
p = plot(xc, T0exp, label="explicit", 
        xlabel="x [m]", ylabel="T [°C]", 
        title="Temperature after $(round(time / day, digits=1)) days
        Δt = $(round(Δt / Δtexp, digits=2))*Δt_{crit}",
        xlim=(0,L),ylim=(0, Tmagma),layout=(1,2))
plot!(p,xc, T0imp,label="implicit",subplot=1)
plot!(p,xc, T0dc,label="def correction",subplot=1)
plot!(p,xc, T0cnv,label="cnv",subplot=1)
plot!(p,xc, Tana, linestyle=:dash, label="analytical",subplot=1)
plot!(p,xc, εexp, xlabel="x [m]", ylabel="ε",
        title="Error",
        label="ε_exp",xlim=(0,L),ylim=(0,2.0),
        subplot=2)        
plot!(p,xc, εimp, label="ε_imp",subplot=2)      
plot!(p,xc, εdc, label="ε_dc",subplot=2)  
plot!(p,xc, εcnv, label="ε_dc",subplot=2)  

if save_fig == 1
    Plots.frame(anim)
else
    display(p)
end
# ----------------------------------------------------------------------- #
# Time loop ------------------------------------------------------------- #
for n=1:nt
    println("Zeitschritt: ",n,", Time: $(round(time/day, digits=1)) [d]")
    # Explicit, Forward Euler ------------------------------------------- #
    T_ex[2:end-1]   =   T0exp
    T1exp       =   ForwardEuler!( T1exp, T_ex, κ, Δt, nc, Δx, BC )
    # Implicit, Backward Euler ------------------------------------------ #
    T1imp,K     =   BackwardEuler!( T0imp, nc, Δx, κ, Δt, BC, K )

    for iter = 1:niter
        # Residual iteration
        ComputeResiduals!(R, T1dc, T_ex, T0dc, ∂T2∂x2, BC, κ, Δx, Δt)
        @printf("||R|| = %1.4e\n", norm(R)/length(R))            
        norm(R)/length(R) < ϵ ? break : nothing
        # Assemble linear system
        K  = AssembleMatrix!(K,BC,nc,κ,Δx,Δt)
        # Solve for temperature correction: Cholesky factorisation
        Kc = cholesky(K.cscmatrix)
        # Solve for temperature correction: Back substitutions
        δT = -(Kc\R[:])
        #δT = -(K\R[:])            
        # Update temperature            
        T = T .+ δT            
    end        
    
    T1cnv,K1,K2     =   CNV!(T0cnv,nc,κ,Δt,Δx,BC,K1,K2)
    
    # Update temperature and time --------------------------------------- #
    T0exp   =   T1exp
    T0imp   =   T1imp
    T0dc    =   T1dc
    T0cnv   =   T1cnv
    
    time    =   time + Δt
    
    # From lecture notes of TWB
    @. Tana     =   Trock + (Tmagma-Trock)/(sqrt(1+4*time*κ/σ^2))*
                        exp(-(xc-xp)^2/(σ^2 + 4*time*κ))
    
    @. εexp     =   abs((Tana-T0exp)/Tana)*100
    @. εimp     =   abs((Tana-T0imp)/Tana)*100
    @. εdc      =   abs((Tana-T0dc)/Tana)*100
    @. εcnv     =   abs((Tana-T0cnv)/Tana)*100
        
    # Plot solution
    if n == 1 || n % 5 == 0 || n == nt
        # Subplot 1 ---
        p = plot(xc, T1exp, label="explicit",
                xlim=(0,L),ylim=(0,1300),
                xlabel="x [m]",ylabel="T [°C]",
                title="Temperature after $(round(time / day, digits=1)) days
        Δt = $(round(Δt / Δtexp, digits=2))*Δt_{crit}",
                layout=(1,2))
        plot!(p, xc, T1imp,linestyle=:dash, label="implicit",subplot=1)
        plot!(p, xc, T1dc,linestyle=:dash, label="def correction",subplot=1)
        plot!(p, xc, T1cnv,linestyle=:dash, label="CNV",subplot=1)
        plot!(p, xc, Tana, linestyle=:dash, label="analytical",subplot=1)    
        # Subplot 2 ---
        plot!(p,xc, εexp, label="ε_exp",
            xlim=(0,L),ylim=(0,2.0),
            xlabel="x [m]",ylabel="ε [%]",
            title="Error",
            subplot=2)
        plot!(p, xc, εimp,linestyle=:dash, label="ε_imp",subplot=2)
        plot!(p, xc, εdc,linestyle=:dash, label="ε_dc",subplot=2)
        plot!(p, xc, εcnv,linestyle=:dash, label="ε_CNV",subplot=2)                
        # Display the plots ---    
        if save_fig == 1
            Plots.frame(anim)
        else
            display(p)
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
foreach(rm, filter(endswith(".png"), readdir(path,join=true)))
# ----------------------------------------------------------------------- #
end
# Call function --------------------------------------------------------- #
Heat_1D_discretization()
# ----------------------------------------------------------------------- #


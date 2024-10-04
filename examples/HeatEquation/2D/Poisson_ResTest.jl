using GeoModBox.HeatEquation.TwoD, ExtendableSparse, Plots, Statistics
using GLM, DataFrames

function Poisson_ResTest()

# Physical Parameters --------------------------------------------------- #
P       = ( 
    L       =   4.0e3,      #   [m]
    H       =   2.0e3,      #   [m]
    k       =   5.6,        #   Conductivity, W/m/K
    # Define the region of the anomaly
    Wcave   =   200.0,      # Width [ m ]
    Hcave   =   200.0,      # Thickness [ m ]
    Dcave   =   1.0e3,      # Depth of center [ m ]
    Xcave   =   2.0e3,      # x-position of center [ m ]
    Q       =   0.3         # volumetric heat production rate [ W/m³ ]; Q = rho*H
)
# ----------------------------------------------------------------------- #
# Initialize grid spacing ----------------------------------------------- #
Δ       = (
    x       =   [0.0],
    y       =   [0.0]
)
# ----------------------------------------------------------------------- #
# Boundary conditions --------------------------------------------------- #
BC      =   (
    type    = (W=:Dirichlet, E=:Dirichlet, N=:Dirichlet, S=:Dirichlet),
    val     = (W=:0.0,E=:0.0,N=:0.0,S=:0.0)
)
# ----------------------------------------------------------------------- #
# Define the numer of differen resolutions ------------------------------ #
# Maximum resolution is defined as nx = n*40, ny = n*20 ---
n       =   20
# Define statistical parameters for visualization ---
ST      =   (
    x       =   zeros(n),       # Reziproc resolution 1/nx/ny
    Tmax    =   zeros(n),       # Maximum temperature
    Tmean   =   zeros(n)        # Mean temperature
)
# ----------------------------------------------------------------------- #
# Loop over the resolutions --------------------------------------------- #
for k = 1:n
    println("k = ",k)
    # Numerical Parameters ---------------------------------------------- #
    NC      = (
        x       =   k*40,        # Gitterpunkte in x-Richtung, Spalten
        y       =   k*20         # Gitterpunkte in y-Richtung, Zeilen    
        
    )       
    
    Δ.x[1]      =   P.L/NC.x
    Δ.y[1]      =   P.H/NC.y
    
    ST.x[k]     =   1 / NC.x / NC.y
    # ------------------------------------------------------------------- #
    # Generate the grid ------------------------------------------------- #
    x       = (
        c       =   LinRange(0.0 + Δ.x[1]/2.0, P.L - Δ.x[1]/2.0, NC.x),
    )
    y       = (
        c       =   LinRange(-P.H + Δ.y[1]/2.0, 0.0 - Δ.y[1]/2.0, NC.y),
    )
    # ------------------------------------------------------------------- #
    # Initialcondition -------------------------------------------------- #
    D       = ( 
        Q       =   zeros(NC...),           # (row,col) 
        T       =   zeros(NC...),
    )
    # Heat production rate in the anomaly ---
    for i = 1:NC.x, j = 1:NC.y
        if x.c[i] >= (P.Xcave-P.Wcave/2.0) && x.c[i] <=(P.Xcave+P.Wcave/2.0) && 
            y.c[j] >= -P.Dcave-P.Hcave/2.0 && y.c[j] <= -P.Dcave+P.Hcave/2.0 
            D.Q[i,j]    = P.Q
        end
    end
    # ------------------------------------------------------------------- #
    # Linear System of Equations ---------------------------------------- #
    Num     =   (T=reshape(1:NC.x*NC.y, NC.x, NC.y),)
    ndof    =   maximum(Num.T)
    K       =   ExtendableSparseMatrix(ndof,ndof)
    rhs     =   zeros(ndof)
    # ------------------------------------------------------------------- #
    # Solve equation ---------------------------------------------------- #
    Poisson!(D,NC,P,BC,Δ,K,rhs,Num)
    # ------------------------------------------------------------------- #
    ST.Tmax[k]      =   maximum(D.T[:])
    ST.Tmean[k]     =   mean(D.T[:])
end
# Linear fit ------------------------------------------------------------ #
df_max      =   DataFrame(x = ST.x, Tmax = ST.Tmax)
df_mean     =   DataFrame(x = ST.x, Tmean = ST.Tmean)
# Fit linear models
linfitMAX   =   lm(@formula(Tmax ~ x), df_max)
linfitMEAN  =   lm(@formula(Tmean ~ x), df_mean)
# Extract coefficients
#coef1       =   coef(linfitMAX)
#coef2       =   coef(linfitMEAN)
# Calculate fitted values
linfit1     =   predict(linfitMAX)
linfit2     =   predict(linfitMEAN)
# ----------------------------------------------------------------------- #
# Plot solution --------------------------------------------------------- #
# Subplot 1 ---
p = scatter(ST.x, ST.Tmean, marker=:circle, markersize=4,label="",
        xlabel="1/nx/ny",ylabel="⟨T⟩",title="Resolution Test",
        xaxis=:log,
        layout=(2,1))
plot!(p,ST.x, linfit2, color="black", label="", 
            linestyle=:dash, linewidth=2)
## Add text for fitted minimum temperature
#annotate!(p,1e-4, 121.8, "T_{fit,mean} = $(coef2[1])")

# subplot 2 ---
scatter!(p,ST.x, ST.Tmax, marker=:circle, markersize=4,label="",
        xlabel="1/nx/ny",ylabel="T_{max}",
        subplot=2)
plot!(p,ST.x, linfit1, color="black", label="", 
            linestyle=:dash, linewidth=2,subplot=2)

# Show the plot
display(p)
savefig("./examples/HeatEquation/2D/Results/Poisson_ResTest.png")
# ----------------------------------------------------------------------- #
end

Poisson_ResTest()

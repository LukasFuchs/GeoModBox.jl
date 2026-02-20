using Plots, GeoModBox.HeatEquation.TwoD, ExtendableSparse
using Statistics, Printf, LinearAlgebra
using TimerOutputs

function Gaussian_Diffusion()
to      =   TimerOutput()
Schema  =   ["explicit","implicit","CNA"]
ns          =   size(Schema,1)
nrnxny      =   6
save_fig    =   1
# Physical Parameters ------------------------------------------------ #
P       = ( 
    L       =   200e3,          #   Length [ m ]
    H       =   200e3,          #   Height [ m ]
    k       =   3,              #   Thermal Conductivity [ W/m/K ]
    cp      =   1000,           #   Specific Heat Capacity [ J/kg/K ]
    ρ       =   3200,           #   Density [ kg/m^3 ]
    K0      =   273.15,         #   Kelvin at 0 °C
    Q0      =   0               #   Heat production rate
)
P1      = (
    κ       =   P.k/P.ρ/P.cp,   #   Thermal Diffusivity [ m^2/s ] 
    Tamp    =   500,            #   Temperaturamplitude [K]
    σ       =   20e3,           #   
    Xc      =   0.0,            #   x-Coordinate of the Anomalycenter
    Zc      =   0.0             #   y-Coordinate of the Anomalycenter
)
P       =   merge(P,P1)
# -------------------------------------------------------------------- #
# Statistical Parameter ---------------------------------------------- #
St      = (
    ε           =   zeros(size(Schema,1),nrnxny),    
    nxny        =   zeros(size(Schema,1),nrnxny),
    Tmax        =   zeros(size(Schema,1),nrnxny),
    Tmean       =   zeros(size(Schema,1),nrnxny),
    Tanamax     =   [0.0],
    Tanamean    =   [0.0]
)
# -------------------------------------------------------------------- #
# Loop over different discretization schemes ------------------------- #
@timeit to "Discretization Loop" begin
for m = 1:ns
    FDSchema = Schema[m]
    display(FDSchema)
    @timeit to "Resolution Loop" begin
    for l = 1:nrnxny
        @timeit to "Ini" begin
        # Numerical Parameters --------------------------------------- #
        NC  = (
            x       =   l*20,       #   Number of Centroids in x
            y       =   l*20        #   Number of Centroids in y
        )
        Δ   = (
            x       =   P.L/NC.x,   #   Grid spacing in x
            y       =   P.H/NC.y    #   Grid Spacing in y
        )
        display(string("nx = ",NC.x,", ny = ",NC.y))
        # ------------------------------------------------------------ #      
        # Grid coordinates ------------------------------------------- #
        x       = (
            c       =   LinRange(-P.L/2+ Δ.x/2.0, P.L/2 - Δ.x/2.0, NC.x),
        )
        y       = (
            c       =   LinRange(-P.H/2 + Δ.y/2.0, P.H/2 - Δ.y/2.0, NC.y),
        )
        # ------------------------------------------------------------ #
        # Time Parameters -------------------------------------------- #
        T       = (
            year        =   365.25*3600*24,     #   Seconds per year
            Δfac        =   1.0,                #   Factor for Explicit Stability Criterion
        )
        T1      = (
            tmax        =   10 * 1e6 * T.year,  #   Maximum Time in [ s ]
            Δ           =   [0.0]            
        )
        T       =   merge(T,T1)
        T.Δ[1]  =   T.Δfac * (1.0 / ( 2.0 * P.κ * ( 1 /Δ.x^2 + 1 / Δ.y^2 )))
        
        nt      =   ceil(Int,T.tmax/T.Δ[1])     #   Number of Time Steps
        time    =   zeros(1,nt)
        # ------------------------------------------------------------ #
        # Initial Conditions  ---------------------------------------- #
        D       = (
            Q           =   zeros(NC...),
            T           =   zeros(NC...),
            T0          =   zeros(NC...),
            T_ex        =   zeros(NC.x+2,NC.y+2),
            T_ex0       =   zeros(NC.x+2,NC.y+2),
            Tana        =   zeros(NC...),
            RMS         =   zeros(1,nt),
            εT          =   zeros(NC...),
            Tmax        =   zeros(1,nt),
            Tmean       =   zeros(1,nt),
            Tmaxa       =   zeros(1,nt),
            Tprofile    =   zeros(NC.y,nt),
            Tprofilea   =   zeros(NC.y,nt),
            ρ           =   zeros(NC...),
            cp          =   zeros(NC...)            
        )
        @. D.ρ  =   P.ρ
        # Initial conditions
        AnalyticalSolution2D!(D.T, x.c, y.c, time[1], (T0=P.Tamp,K=P.κ,σ=P.σ))
        @. D.Tana                   =   D.T
        @. D.T0                     =   D.T
        D.T_ex[2:end-1,2:end-1]     .=  D.T
    
        D.Tprofile[:,1]     .=  (D.T[convert(Int,NC.x/2),:] + 
                                    D.T[convert(Int,NC.x/2)+1,:]) / 2
        D.Tprofilea[:,1]    .=  (D.Tana[convert(Int,NC.x/2),:] + 
                                    D.Tana[convert(Int,NC.x/2)+1,:]) / 2
        # Heat production rate ---
        @. D.Q          = P.Q0
        # Boundary Conditions ---------------------------------------- #
        BC     = (type    = (W=:Dirichlet, E=:Dirichlet, 
                                N=:Dirichlet, S=:Dirichlet),
                    val     = (W=D.Tana[1,:],E=D.Tana[end,:],
                                N=D.Tana[:,end],S=D.Tana[:,1]))
        # ------------------------------------------------------------ #
        niter       =   10
        ϵ           =   1e-20
        @. D.ρ      =   P.ρ
        @. D.cp     =   P.cp
        k           =   (x=zeros(NC.x+1,NC.x), y=zeros(NC.x,NC.x+1))
        @. k.x      =   P.k
        @. k.y      =   P.k
        Num         =   (T=reshape(1:NC.x*NC.y, NC.x, NC.y),)
        ndof        =   maximum(Num.T)
        K           =   ExtendableSparseMatrix(ndof,ndof)
        R           =   zeros(NC...)
        ∂T          =   (∂x=zeros(NC.x+1, NC.x), ∂y=zeros(NC.x, NC.x+1))
        q           =   (x=zeros(NC.x+1, NC.x), y=zeros(NC.x, NC.x+1),
                            x0=zeros(NC.x+1, NC.x), y0=zeros(NC.x, NC.x+1),)
        if FDSchema == "implicit"
            C = 0
        elseif FDSchema == "CNA"
            C = 0.5
        elseif FDSchema == "explicit"
            C = 1
        end
        end
        @timeit to "Time Loop" begin
        # Time Loop -------------------------------------------------- #
        for n = 1:nt
            if n>1
                @timeit to "Solution" begin
                for iter = 1:niter
                    # Evaluate residual
                    ComputeResiduals2D!(R, D.T, D.T_ex, D.T0, D.T_ex0, D.Q, ∂T, 
                            q, D.ρ, D.cp, k, BC, Δ, T.Δ[1];C)
                    # @printf("||R|| = %1.4e\n", norm(R)/length(R))
                    norm(R)/length(R) < ϵ ? break : nothing
                    # Assemble linear system
                    K  = AssembleMatrix2D(D.ρ, D.cp, k, BC, Num, NC, Δ, T.Δ[1];C)
                    # Solve for temperature correction: Cholesky factorisation
                    Kc = cholesky(K.cscmatrix)
                    # Solve for temperature correction: Back substitutions
                    δT = -(Kc\R[:])
                    # Update temperature
                    @. D.T += δT[Num.T]
                end
                D.T0    .= D.T
                end
                time[n]     =   time[n-1] + T.Δ[1]
                if time[n] > T.tmax 
                    T.Δ[1]  =   T.tmax - time[n-1]
                    time[n] =   time[n-1] + T.Δ[1]
                end                
                # Exact solution on cell centroids
                AnalyticalSolution2D!(D.Tana, x.c, y.c, time[n], (T0=P.Tamp,K=P.κ,σ=P.σ))
                # Exact solution on cell boundaries
                BoundaryConditions2D!(BC, x.c, y.c, time[n], (T0=P.Tamp,K=P.κ,σ=P.σ)) 
            end
            # Maximum and Mean Temperature with time ---
            D.Tmax[n]   =   maximum(D.T)
            D.Tmean[n]  =   mean(D.T)
            # Vertical Profile along the Center of the Domain ---
            D.Tprofile[:,n]     .=  (D.T[convert(Int,NC.x/2),:] + 
                                        D.T[convert(Int,NC.x/2)+1,:]) / 2
            D.Tprofilea[:,n]    .=  (D.Tana[convert(Int,NC.x/2),:] + 
                                        D.Tana[convert(Int,NC.x/2)+1,:]) / 2
            # Deviation from the Analytical Solution ---
            @. D.εT     =   (D.Tana - D.T)
            # RMS ---
            D.RMS[n]    =   sqrt(sum(D.εT.^2)/(NC.x*NC.y))
        end        
        display("Time loop finished ...")
        display("-> Use new grid size...")
        end
        # ------------------------------------------------------------ #
        # Statistical Values for Each Scheme and Resolution ---
        St.ε[m,l]       =   maximum(D.RMS[:])
        St.nxny[m,l]    =   1/NC.x/NC.y
        St.Tmax[m,l]    =   D.Tmax[nt]
        St.Tmean[m,l]   =   D.Tmean[nt]
        St.Tanamax[1]   =   maximum(D.Tana)
        St.Tmean[1]     =   mean(D.Tana)
        # ------------------------------------------------------------ #
    end
    end
end
end
# Visualize Statistical Values --------------------------------------- #
q   =   plot(0,0,layout=(1,3))
for m = 1:ns
    plot!(q,St.nxny[m,:],St.ε[m,:],
                marker=:circle,markersize=3,label=Schema[m],
                xaxis=:log,yaxis=:log,
                xlabel="1/nx/ny",ylabel="ε_{T}",layout=(1,3),
                subplot=1)
    plot!(q,St.nxny[m,:],St.Tmax[m,:],
                marker=:circle,markersize=3,label="",
                xaxis=:log,
                xlabel="1/nx/ny",ylabel="T_{max}",
                subplot=2)
    plot!(q,St.nxny[m,:],St.Tmean[m,:],
                marker=:circle,markersize=3,label="",
                xaxis=:log,
                xlabel="1/nx/ny",ylabel="⟨T⟩",
                subplot=3)
    display(q)
end
# --------------------------------------------------------------------- #
# Save Final Figure --------------------------------------------------- #
if save_fig == 1
    savefig(q,"./examples/DiffusionEquation/2D/Results/Gaussian_ResTest_General_variable_k.png")
end
# --------------------------------------------------------------------- #
display(to)
end

Gaussian_Diffusion()


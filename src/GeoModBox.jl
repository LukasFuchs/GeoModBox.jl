module GeoModBox

    using Statistics: mean
    using ExtendableSparse, LinearAlgebra

    module HeatEquation        

        module TwoD
            # Handle analytical solution for 2D Diffusion2D_Gaussian
            include("./HeatEquation/AnalyticsDiffusion2D.jl")
            export AnalyticalSolution!, BoundaryConditions!	    

            # Implicit solver
            include("./HeatEquation/BackwardEuler.jl")
            export ComputeResiduals!, AssembleMatrix 

            # Poisson solver
            include("./HeatEquation/PoissonSolvers.jl")
            export Poisson! 
        end

        module OneD
            # 1D solver
            include("./HeatEquation/1Dsolvers.jl")
            export ForwardEuler!, BackwardEuler!,
                ComputeResiduals!, AssembleMatrix!, CNA!
        end        
    end
end 

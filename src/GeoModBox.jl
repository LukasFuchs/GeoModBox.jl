module GeoModBox

    using Statistics: mean
    using ExtendableSparse, LinearAlgebra

    include("./Structures.jl")
    export Geometry, Physics, 
            GridSpacing,
            DataFields, TimeParameter

    module Scaling
   
        include("./Scaling.jl")
        export Constants, ScalingConstants!, ScaleParameters!

    end

    module HeatEquation            

        module TwoD
            # Handle analytical solution for 2D Diffusion2D_Gaussian ---
            include("./HeatEquation/AnalyticsDiffusion2D.jl")
            export AnalyticalSolution2D!, BoundaryConditions2D!	    

            # 2D solver ---
            include("./HeatEquation/2Dsolvers.jl")
            export ComputeResiduals2D!, AssembleMatrix2D, BackwardEuler2Dc!, 
                ForwardEuler2Dc!, CNA2Dc!, ADI2Dc!,
                Poisson2Dc!, Poisson2D!
        end

        module OneD
            # 1D solver ---
            include("./HeatEquation/1Dsolvers.jl")
            export ForwardEuler1Dc!, BackwardEuler1Dc!,
                ComputeResiduals1Dc!, AssembleMatrix1Dc!, CNA1Dc!,
                ForwardEuler1D!
        end                
    end

    module AdvectionEquation
    
        module OneD
            # 1D solver ---
            include("./AdvectionEquation/1Dsolvers.jl")
            export upwind1D!, lax1D!, slf1D!, semilag1D!, RK4O1D!            
        end

        module TwoD
            # 2D solver --- 
            include("./AdvectionEquation/2Dsolvers.jl")
            export upwindc2D!, slfc2D!, semilagc2D!
        end
    end

    module InitialCondition

        # Initial Conditions ---
        include("./InitialCondition/2Dini.jl")        
        export IniVelocity!, IniTemperature!, IniPhase!
    end

    module Tracers

        module OneD
            include("./Tracers/1Dsolvers.jl")
            export Itp1D_Centers2Markers!, Itp1D_Markers2Centers!
        end

        module TwoD
            include("./Tracers/2Dsolvers.jl")
            export TMarkers, Markers, IniTracer2D, 
                    VxFromVxNodes, VyFromVyNodes, VxVyFromPrNodes,
                    FromCtoM, CountMPC, Markers2Cells, AdvectTracer2D,
                    Markers2Vertices
        end
    end

    module MomentumEquation
        
        module OneD
            include("./MomentumEquation/1Dsolvers.jl")
            export Stokes_1D_direct, ComputeStokesResiduals1D!, 
                    AssembleStokesMatrix1D
        end

        module TwoD
            include("./MomentumEquation/2Dsolvers.jl")
            export Assemblyc, updaterhsc, Residuals2Dc!, 
                    Assembly, updaterhs, Residuals2D!
        end
    end

end 

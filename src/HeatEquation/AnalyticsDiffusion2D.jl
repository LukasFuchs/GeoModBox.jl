using ExactFieldSolutions, StaticArrays

@doc raw"""
    AnalyticalSolution!(Te, x, y, t)

Calls 2D analytical solution for a 2D diffusion problem:

    Te     : 2D matrix 
    x      : x coordinate array
    y      : y coordinate array
    t      : time
# Examples
```julia-repl
julia> a = zeros(2,2); x = [0, 1]; y = [0, 1];
julia> AnalyticalSolution!(a, x, y, 0.0)
2×2 Matrix{Float64}:
 100.0          6.04202e-67
   6.04202e-67  3.6506e-135
```
"""
@views function AnalyticalSolution!(Te, x, y, t)
    for i in axes(Te,1), j in axes(Te,2)
        X = @SVector([x[i], y[j], t])
        s = Diffusion2D_Gaussian( X )
        Te[i,j] = s.u
    end
end

@views function BoundaryConditions!(BC, x, y, t) 
    # Box boundaries coordinates
    xmin = x[1]   - (x[2]-x[1])/2
    xmax = x[end] + (x[2]-x[1])/2
    ymin = y[1]   - (y[2]-y[1])/2
    ymax = y[end] + (y[2]-y[1])/2
    # Loop over W/E sides
    for j in axes(BC.val.W,1)
        X = @SVector([xmin, y[j], t])
        s = Diffusion2D_Gaussian( X )
        BC.val.W[j] = s.u
        X = @SVector([xmax, y[j], t])
        s = Diffusion2D_Gaussian( X )
        BC.val.E[j] = s.u
    end
    # Loop over S/N sides
    for i in axes(BC.val.S,1)
        X = @SVector([x[i], ymin, t])
        s = Diffusion2D_Gaussian( X )
        BC.val.S[i] = s.u
        X = @SVector([x[i], ymax, t])
        s = Diffusion2D_Gaussian( X )
        BC.val.N[i] = s.u
    end
end
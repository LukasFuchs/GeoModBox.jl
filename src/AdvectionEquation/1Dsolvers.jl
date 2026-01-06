using Dierckx

"""
    RK4O1D!( x, Δt, vx, xmin, xmax )    

Function to advect a tracer in 1-D for one time step using fourth order Runge-Kutta. 

The final position of the tracer is defined by: 

    x   .+  =   (x1 + 2.0 * (x2 + x3) + x4) / 6.0 

where 

    x1  =   Δt * vx 
    x2  =   Δt * vx
    x3  =   Δt * vx
    x4  =   Δt * vx 

Examples
========
```julia
julia> xminx,xmax = 0,100
(0, 100)

julia> x = collect(xmin:0.1:xmax)
1001-element Vector{Float64}:
   0.0
   0.1
   0.2
   0.3
   0.4
   ⋮
  99.7
  99.8
  99.9
 100.0

 julia> x = RK4O1D!(x,0.1,1,xmin,xmax)
```

"""
function RK4O1D!( x, Δt, vx, xmin, xmax )    

    x1  =   Δt * vx 
    x2  =   Δt * vx
    x3  =   Δt * vx
    x4  =   Δt * vx 

    x    .+=   (x1 + 2.0 * (x2 + x3) + x4) / 6.0 

    x[x.>xmax] .= xmin .+ abs.(x[x.>xmax] .- xmax)
    x[x.<xmin] .= xmax .- abs.(x[x.<xmin] .- xmin)
    
end

"""
    upwind1D!( A, vx, Δt, Δx )

Function to advect a tracer using the upwind method in one dimension. 

    A   : 1-D Array of the advected property
    vx  : Horizontal velocity
    Δt  : Time step 
    Δx  : Grid spacing

"""
function upwind1D!( A, vx, Δt, Δx )
    Aold    =   zeros(size(A))

    Aold    .=  A

    if vx > 0
        @. A[2:end-1] = 
            Aold[2:end-1] - vx*Δt/Δx*( Aold[2:end-1] - Aold[1:end-2] )
    elseif vx < 0
        @. A[2:end-1] = 
            Aold[2:end-1] - vx*Δt/Δx*( Aold[3:end] - T[2:end-1] ) 
    end
end

"""
    lax1D!( A, vx, Δt, Δx )

Function to advect a property using the Lax-Friedrich method in one dimension. 

    A   : 1-D array of the advected property
    vx  : Horizontal velocity
    Δt  : Time step
    Δx  : Grid spacing
"""
function lax1D!( A, vx, Δt, Δx )
    Aold    =   zeros(size(A))

    Aold    .=  A

    @. A[2:end-1] = ( Aold[3:end] + Told[1:end-2])/2 -
            (vx*Δt/2/Δx) * ( Aold[3:end] - Aold[1:end-2])
end

"""
    slf1D!( A, Aold2, vx, Δt, Δx )

Function to advect a property using the staggered leaped frog method in one dimension. 

    A       : 1-D array of the advected property
    Aold2   : 1-D array of the advected property at the previous time step
    vx      : Horizontal velocity
    Δt      : Time step
    Δx      : Grid spacing
"""
function slf1D!( A, Aold2, vx, Δt, Δx )
    # Aold2 is old time step

    Aold    =   zeros(size(A))

    Aold    .=  A

    @. T[2:end-1] = 
        Aold2[2:end-1] - vx*Δt/Δx * ( Aold[3:end] - Aold[1:end-2])

end

"""
    semilag1D!( A, xc, vx, Δt, Δx )

Function to advect a property using the semi-lagrangian method in one dimension. 

    A   : 1-D array of the advected property
    xc  : x-coordinates of each centroid
    vx  : Horizontal velocity
    Δt  : Time step
    Δx  : Horizontal velocity
"""
function semilag1D!( A, xc, vx, Δt, Δx )
    X       =   zeros(size(xc))
    Asl     =   zeros(size(A))

    @. X    =   xc - Δt*vx
        
    spl     =   Spline1D(xc, A; k=1)
    Asl     =   spl.(X)
    A       .=  Asl     
end
using Dierckx

"""
    RK4O1D!( x, Δt, vx, xmin, xmax )    

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
"""
function lax1D!( A, vx, Δt, Δx )
    Aold    =   zeros(size(A))

    Aold    .=  A

    @. A[2:end-1] = ( Aold[3:end] + Told[1:end-2])/2 -
            (vx*Δt/2/Δx) * ( Aold[3:end] - Aold[1:end-2])
end

"""
    slf1D!( A, Aold2, vx, Δt, Δx )
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
"""
function semilag1D!( A, xc, vx, Δt, Δx )
    X       =   zeros(size(xc))
    Asl     =   zeros(size(A))

    @. X    =   xc - Δt*vx
        
    spl     =   Spline1D(xc, A; k=1)
    Asl     =   spl.(X)
    A       .=  Asl     
end
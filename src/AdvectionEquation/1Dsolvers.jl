using Dierckx

"""
    RK4O1D!(x, Δt, vx, xmin, xmax)

Advects Lagrangian tracers in one dimension using the classical fourth-order
Runge–Kutta (RK4) method.

The tracer trajectories are integrated explicitly over one time step assuming
a constant advection velocity. Under this assumption, all four Runge–Kutta
stages evaluate the same velocity, reducing the method to an exact integration
of the particle trajectories for uniform flow. Tracers leaving one side of the
computational domain are reinserted at the opposite side, resulting in periodic
boundary conditions.

# Arguments

    x        : Vector containing the tracer coordinates.
    Δt       : Time step.
    vx       : Constant horizontal advection velocity.
    xmin     : Minimum x-coordinate of the computational domain.
    xmax     : Maximum x-coordinate of the computational domain.

# Notes

The routine assumes a spatially and temporally constant velocity field during
the current time step. Periodic boundary conditions are enforced by wrapping
tracers that leave the computational domain to the opposite boundary.

"""
function RK4O1D!( x, Δt, vx, xmin, xmax )    

    x1  =   Δt * vx 
    x2  =   Δt * vx
    x3  =   Δt * vx
    x4  =   Δt * vx 

    x    .+=   (x1 + 2.0 * (x2 + x3) + x4) / 6.0 

    L = xmax - xmin
    @. x = xmin + mod(x - xmin, L)

    # x[x.>xmax] .= xmin .+ abs.(x[x.>xmax] .- xmax)
    # x[x.<xmin] .= xmax .- abs.(x[x.<xmin] .- xmin)
    
end

"""

    upwind1D!(A, vx, Δt, Δx)

Advects a scalar property using a first-order upwind finite-difference scheme
on a one-dimensional grid.

The spatial derivative is approximated using a one-sided finite difference,
with the upwind direction selected according to the sign of the constant
advection velocity. The method is explicit, robust, and monotonic, although
it introduces numerical diffusion that increases with the Courant number and
the number of time steps.

# Arguments

    A        : Vector containing the advected property.
    vx       : Constant horizontal advection velocity.
    Δt       : Time step.
    Δx       : Grid spacing.

# Notes

The advection update is performed explicitly using first-order upwind finite
differences. The implementation assumes a spatially constant velocity during
the current time step.

For stability, the time step should satisfy the
Courant–Friedrichs–Lewy (CFL) condition.

"""
function upwind1D!( A, vx, Δt, Δx )

    Aold    =  copy(A)

    if vx > 0
        @. A[2:end-1] = 
            Aold[2:end-1] - vx*Δt/Δx*( Aold[2:end-1] - Aold[1:end-2] )
    elseif vx < 0
        @. A[2:end-1] = 
            Aold[2:end-1] - vx*Δt/Δx*( Aold[3:end] - Aold[2:end-1] ) 
    end
end

"""
    lax1D!(A, vx, Δt, Δx)

Advects a scalar property using the first-order Lax–Friedrichs finite-difference
scheme on a one-dimensional grid.

The Lax–Friedrichs method replaces the property at each grid point by the
average of its neighboring values before applying a central-difference
approximation of the advection term. This additional averaging introduces
artificial numerical diffusion, which stabilizes the explicit scheme and
suppresses non-physical oscillations.

# Arguments

    A        : Vector containing the advected property.
    vx       : Constant horizontal advection velocity.
    Δt       : Time step.
    Δx       : Grid spacing.

# Notes

The advection update is performed explicitly using a combination of neighbor
averaging and central finite differences. The implementation assumes a
spatially constant velocity during the current time step.

Compared to the first-order upwind scheme, the Lax–Friedrichs method is
symmetric with respect to the flow direction but generally introduces stronger
numerical diffusion.

For stability, the time step should satisfy the
Courant–Friedrichs–Lewy (CFL) condition.
"""
function lax1D!( A, vx, Δt, Δx )

    Aold    =  copy(A)

    @. A[2:end-1] = ( Aold[3:end] + Aold[1:end-2])/2 -
            (vx*Δt/2/Δx) * ( Aold[3:end] - Aold[1:end-2])
end

"""
    slf1D!(A, Aold2, vx, Δt, Δx)

Advects a scalar property using the staggered leapfrog (SLF) finite-difference
scheme on a one-dimensional grid.

The staggered leapfrog method is an explicit second-order accurate advection
scheme in both space and time. The spatial derivative is approximated using a
central finite difference, while the leapfrog time integration combines the
solution from the current and previous time steps. Compared to first-order
schemes, the staggered leapfrog method exhibits significantly lower numerical
diffusion but may develop dispersive oscillations near sharp gradients.

# Arguments

    A        : Vector containing the advected property at the current time step.
    Aold2    : Vector containing the advected property from the previous time step.
    vx       : Constant horizontal advection velocity.
    Δt       : Time step.
    Δx       : Grid spacing.

# Notes

The advection update is performed explicitly using a central finite-difference
approximation of the spatial derivative and leapfrog time integration. The
implementation assumes a spatially constant velocity during the current time
step.

For stability, the time step should satisfy the
Courant–Friedrichs–Lewy (CFL) condition.
"""
function slf1D!( A, Aold2, vx, Δt, Δx )
    Aold    =  copy(A)

    @. A[2:end-1] = 
        Aold2[2:end-1] - vx*Δt/Δx * ( Aold[3:end] - Aold[1:end-2])

    Aold2 .= Aold

end

"""
    semilag1D!( A, xc, vx, Δt )

Advects a scalar property using a semi-Lagrangian method on a one-dimensional
grid.

The departure point of each centroid is calculated by tracing the trajectory
backward over one time step using the prescribed horizontal velocity. The
advected property is then evaluated at the departure points using linear spline
interpolation.

# Arguments

    A        : Vector containing the advected property.
    xc       : Vector containing the centroid x-coordinates.
    vx       : Horizontal advection velocity.
    Δt       : Time step.
    Δx       : Grid spacing.

# Notes

The implementation assumes that the velocity remains constant during the
current time step. Linear interpolation is used to evaluate the property at the
departure points.

Unlike explicit Eulerian advection schemes, the semi-Lagrangian method is not
restricted by the conventional advective CFL stability condition. However,
the time step still affects the trajectory and interpolation accuracy.
"""
function semilag1D!( A, xc, vx, Δt )
    X = similar(xc)
    @. X = xc - Δt * vx

    spl = Spline1D(xc, A; k=1)
    A  .= spl.(X)

    return nothing
end
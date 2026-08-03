# -------------------------------------------------------------------- #
# 2D solver for the advection equation 
# -------------------------------------------------------------------- #
using Interpolations

"""
    upwindc2D!(P, P_ex, vxc, vyc, NC, Δt, Δx, Δy)

Advects a scalar property using a first-order upwind finite-difference scheme
on a two-dimensional finite-difference grid.

The spatial derivatives are evaluated using one-sided finite differences
selected according to the sign of the local velocity components. This
upwind discretization provides a robust and monotonic solution, although
it introduces numerical diffusion, particularly when the flow is not aligned
with the computational grid.

# Arguments

    P        : 2-D array containing the property defined on the centroids.
    P_ex     : 2-D array containing the property including ghost nodes.
    vxc      : Horizontal velocity defined on the centroids.
    vyc      : Vertical velocity defined on the centroids.
    NC       : Structure or tuple containing the number of centroids.
    Δt       : Time step.
    Δx       : Horizontal grid spacing.
    Δy       : Vertical grid spacing.

# Notes

The advection update is performed explicitly using first-order upwind
differences in both coordinate directions. After each advection step,
the interior values are copied back to the extended array containing the
ghost nodes, allowing the boundary conditions to be updated before the
next time step.

For stability, the time step should satisfy the Courant–Friedrichs–Lewy
(CFL) condition.
"""
function upwindc2D!(P,P_ex,vxc,vyc,NC,Δt,Δx,Δy)

    @. P     =  P_ex[2:(NC.x+1),2:(NC.y+1)] - 
            (vxc>0)*(vxc*Δt/Δx*(P_ex[2:(NC.x+1),2:(NC.y+1)] - P_ex[1:NC.x,2:(NC.y+1)])) - 
            (vxc<0)*(vxc*Δt/Δx*(P_ex[3:(NC.x+2),2:(NC.y+1)] - P_ex[2:(NC.x+1),2:(NC.y+1)])) - 
            (vyc>0)*(vyc*Δt/Δy*(P_ex[2:(NC.x+1),2:(NC.y+1)] - P_ex[2:(NC.x+1),1:NC.y])) - 
            (vyc<0)*(vyc*Δt/Δy*(P_ex[2:(NC.x+1),3:(NC.y+2)] - P_ex[2:(NC.x+1),2:(NC.y+1)]))
    # Update extended temperature field ------------------------------ #
    P_ex[2:(NC.x+1),2:(NC.y+1)]     .=  P
    # ---------------------------------------------------------------- #
end

"""
    slfc2D!(P, P_ex, P_exo, vxc, vyc, NC, Δt, Δx, Δy)

Advects a scalar property using the staggered leapfrog (SLF) finite-difference
scheme on a two-dimensional finite-difference grid.

The staggered leapfrog method is an explicit second-order accurate advection
scheme in both space and time. Spatial derivatives are approximated using
central finite differences, while the leapfrog time integration combines the
property fields from the current and previous time steps. Compared to the
first-order upwind scheme, the staggered leapfrog method exhibits significantly
lower numerical diffusion but may develop dispersive oscillations in regions
with strong gradients.

# Arguments

    P        : 2-D array containing the property defined on the centroids.
    P_ex     : 2-D array containing the property including ghost nodes at the current time step.
    P_exo    : 2-D array containing the property including ghost nodes from the previous time step.
    vxc      : Horizontal velocity defined on the centroids.
    vyc      : Vertical velocity defined on the centroids.
    NC       : Structure or tuple containing the number of centroids.
    Δt       : Time step.
    Δx       : Horizontal grid spacing.
    Δy       : Vertical grid spacing.

# Notes

The advection update is performed explicitly using central finite differences
for the spatial derivatives and a leapfrog time integration. After each
advection step, the current property field is copied to `P_exo`, and the newly
computed solution is written back to the interior of `P_ex`. The ghost nodes
can subsequently be updated by the appropriate boundary-condition routine
before the next time step.

For stability, the time step should satisfy the Courant–Friedrichs–Lewy
(CFL) condition.
"""
function slfc2D!(P,P_ex,P_exo,vxc,vyc,NC,Δt,Δx,Δy)

    indx    =   2:(NC.x+1)
    indy    =   2:(NC.y+1)

    @. P  =   P_exo[indx,indy] - 
        vxc*Δt/Δx*(P_ex[indx+1,indy] - P_ex[indx.-1,indy]) - 
        vyc*Δt/Δy*(P_ex[indx,indy+1] - P_ex[indx,indy-1])
    # Update extended temperature field ------------------------------ #
    @. P_exo            =  P_ex
    P_ex[indx,indy]     .=  P
    # ---------------------------------------------------------------- #
end

"""
    semilagc2D!(P, P_ex, vxc, vyc, vxo, vyo, x, y, Δt;
                maxiter=10, tol=1e-10, verbose=false)

Advects a scalar property using a second-order semi-Lagrangian method on a
two-dimensional finite-difference grid.

The departure point of each centroid is computed by integrating the particle
trajectory backward in time using an implicit midpoint scheme. The midpoint
location is obtained through a fixed-point iteration using a time-centered
velocity field, defined as the average of the velocity fields from the current
and previous time steps. The velocity field is interpolated linearly during the
iteration, while the transported property is interpolated using cubic spline
interpolation.

# Arguments

    P        : 2-D array containing the property defined on the centroids.
    P_ex     : 2-D array containing the property including ghost nodes.
    vxc      : Horizontal centroid velocity at the current time step.
    vyc      : Vertical centroid velocity at the current time step.
    vxo      : Horizontal centroid velocity at the previous time step.
    vyo      : Vertical centroid velocity at the previous time step.
    x        : Structure or tuple containing the centroid x-coordinates.
    y        : Structure or tuple containing the centroid y-coordinates.
    Δt       : Time step.

# Keyword Arguments

    maxiter  : Maximum number of midpoint iterations (default: 10).
    tol       : Convergence tolerance for the midpoint iteration (default: 1e-10).
    verbose   : Print convergence information if `true` (default: false).

# Notes

The implicit midpoint iteration typically converges within only a few
iterations. Once the departure points have been determined, the transported
property is interpolated from the previous time step onto the departure
locations using cubic spline interpolation.
"""

function semilagc2D!(
        P, P_ex,
        vxc, vyc,
        vxo, vyo,
        x, y, Δt;
        maxiter = 10,
        tol = 1e-10,
        verbose = false,
    )

    # Time-centred velocity
    if isempty(vxo) || isempty(vyo)
        vxcm = copy(vxc)
        vycm = copy(vyc)
    else
        vxcm = @. 0.5 * (vxc + vxo)
        vycm = @. 0.5 * (vyc + vyo)
    end

    # Velocity interpolators
    itp_vx = linear_interpolation(
        (x.c, y.c),
        vxcm;
        extrapolation_bc = Line(),
    )

    itp_vy = linear_interpolation(
        (x.c, y.c),
        vycm;
        extrapolation_bc = Line(),
    )

    # Initial midpoint-velocity estimate
    vxi = copy(vxcm)
    vyi = copy(vycm)

    vxnew = similar(vxi)
    vynew = similar(vyi)

    xp = similar(vxi)
    yp = similar(vyi)

    converged = false

    # Implicit midpoint iteration
    for k = 1:maxiter

        @. xp = x.c2d - 0.5 * Δt * vxi
        @. yp = y.c2d - 0.5 * Δt * vyi

        if !all(isfinite, xp) || !all(isfinite, yp)
            error("Non-finite midpoint coordinates during iteration $k.")
        end

        vxnew .= itp_vx.(xp, yp)
        vynew .= itp_vy.(xp, yp)

        velocity_scale = max(
            maximum(abs.(vxnew)),
            maximum(abs.(vynew)),
            eps(Float64),
        )

        err = max(
            maximum(abs.(vxnew .- vxi)),
            maximum(abs.(vynew .- vyi)),
        ) / velocity_scale

        if !isfinite(err)
            error("Semi-Lagrangian midpoint iteration diverged.")
        end

        vxi .= vxnew
        vyi .= vynew

        if err < tol
            converged = true
            verbose &&
                println("Midpoint iteration converged after $k iterations.")
            break
        end
    end

    if !converged && verbose
        @warn "Midpoint iteration reached maxiter = $maxiter."
    end

    # Full departure points
    @. xp = x.c2d - Δt * vxi
    @. yp = y.c2d - Δt * vyi

    if !all(isfinite, xp) || !all(isfinite, yp)
        error("Semi-Lagrangian departure coordinates are not finite.")
    end

    if minimum(xp) < first(x.ce) ||
       maximum(xp) > last(x.ce)  ||
       minimum(yp) < first(y.ce) ||
       maximum(yp) > last(y.ce)

        error(
            "Semi-Lagrangian departure points lie outside the extended " *
            "property domain. Check the time step, velocity boundary " *
            "conditions, and scalar boundary treatment."
        )
    end

    # Interpolate the transported property
    itp_P = cubic_spline_interpolation(
        (x.ce, y.ce),
        P_ex,
    )

    P .= itp_P.(xp, yp)

    P_ex[2:end-1, 2:end-1] .= P

    return nothing
end

# function semilagc2D!(
#         P, P_ex,
#         vxc, vyc,
#         vxo, vyo,
#         x, y, Δt;
#         maxiter = 10,
#         tol = 1e-10,
#         verbose = false)

#     #--------------------------------------------------------------------------
#     # Time centred velocity
#     #--------------------------------------------------------------------------
#     vxcm = 0.5 .* (vxc .+ vxo)
#     vycm = 0.5 .* (vyc .+ vyo)

#     #--------------------------------------------------------------------------
#     # Interpolators (build ONCE)
#     #--------------------------------------------------------------------------
#     itp_vx = linear_interpolation(
#         (x.c, y.c),
#         vxcm;
#         extrapolation_bc = Line()
#     )

#     itp_vy = linear_interpolation(
#         (x.c, y.c),
#         vycm;
#         extrapolation_bc = Line()
#     )

#     # Cubic interpolation of transported quantity
#     itp_P = cubic_spline_interpolation(
#         (x.ce, y.ce),
#         P_ex
#     )

#     #--------------------------------------------------------------------------
#     # Initial guess (arrival velocity)
#     #--------------------------------------------------------------------------
#     vxi = copy(vxcm)
#     vyi = copy(vycm)

#     xp = similar(vxi)
#     yp = similar(vyi)

#     converged = false

#     #--------------------------------------------------------------------------
#     # Implicit midpoint iteration
#     #--------------------------------------------------------------------------
#     for k = 1:maxiter

#         @. xp = x.c2d - 0.5*Δt*vxi
#         @. yp = y.c2d - 0.5*Δt*vyi

#         #--------------------------------------------------------------
#         # Diagnostics
#         #--------------------------------------------------------------
#         if !all(isfinite, xp) || !all(isfinite, yp)
#             error("Non-finite midpoint coordinates during iteration $k")
#         end

#         # Uncomment for debugging
#         # @assert minimum(xp) > first(x.ce)-1e-6
#         # @assert maximum(xp) < last(x.ce)+1e-6

#         vxnew = itp_vx.(xp, yp)
#         vynew = itp_vy.(xp, yp)

#         err = max(
#             maximum(abs.(vxnew .- vxi)),
#             maximum(abs.(vynew .- vyi))
#         )

#         if !isfinite(err)
#             error("Midpoint iteration diverged.")
#         end

#         vxi .= vxnew
#         vyi .= vynew

#         if err < tol
#             converged = true
#             verbose && println("Midpoint converged after $k iterations.")
#             break
#         end

#     end

#     if !converged && verbose
#         @warn "Midpoint iteration reached maximum iterations."
#     end

#     #--------------------------------------------------------------------------
#     # Full departure point
#     #--------------------------------------------------------------------------
#     @. xp = x.c2d - Δt*vxi
#     @. yp = y.c2d - Δt*vyi

#     if !all(isfinite, xp) || !all(isfinite, yp)
#         error("Departure coordinates are not finite.")
#     end

#     #--------------------------------------------------------------------------
#     # Interpolate scalar
#     #--------------------------------------------------------------------------
#     P .= itp_P.(xp, yp)

#     P_ex[2:end-1,2:end-1] .= P

#     return

# end
# Advection Equation

In most geodynamical applications, material is in motion, and physical
properties such as temperature, density, composition, or accumulated strain
are transported by the flow field. This process is referred to as
*advection*. Although the following discussion focuses primarily on
temperature, the numerical methods described here can also be applied to other
scalar properties.

Depending on the chosen reference frame, advection can be expressed in either
an Eulerian or a Lagrangian formulation.

## Eulerian and Lagrangian Formulations

### Eulerian Formulation

In a fixed Eulerian reference frame, the local rate of change of temperature
is described by the advection equation

$\begin{equation}
\frac{\partial T}{\partial t} = - v_i \frac{\partial T}{\partial x_i},
\end{equation}$

where $T$ is the temperature field and $v_i$ denotes the velocity component
in the $i$-direction. Using vector notation, this equation can equivalently be
written as

$\begin{equation} 
\frac{\partial T}{\partial t} + \mathbf{v}\cdot\nabla T = 0.
\end{equation}$

For incompressible flow,

$\begin{equation}
\nabla\cdot\mathbf{v}=0,
\end{equation}$

the advective and conservative forms of the advection equation are equivalent:

$\begin{equation}
\frac{\partial T}{\partial t} + \nabla\cdot\left(\mathbf{v}T\right) = 0.
\end{equation}$

The advection equation is a first-order hyperbolic partial differential
equation. Information is transported along characteristic trajectories
defined by the velocity field.

Eulerian formulations are solved on a fixed numerical grid using
finite-difference schemes designed to balance accuracy, stability, and
computational efficiency. `GeoModBox.jl` provides several such schemes for
one- and two-dimensional problems.

### Lagrangian Formulation

In a Lagrangian reference frame, the numerical description follows individual
material particles. The evolution of temperature along a particle trajectory
is described by the material derivative

$\begin{equation}
\frac{D T}{D t} = \frac{\partial T}{\partial t} + v_i\frac{\partial T}{\partial x_i}.
\end{equation}$

For pure advection without diffusion or internal heat production,

$\begin{equation}
\frac{D T}{D t}=0,
\end{equation}$

meaning that the temperature carried by an individual material particle
remains constant along its trajectory.

The particle trajectories are obtained by solving the system of ordinary
differential equations

$\begin{equation}
\frac{D x_i}{D t}=v_i,
\end{equation}$

where $x_i$ denotes the particle coordinate in the $i$-direction.

## Discretization Schemes

Although the governing equation is simple in form, its numerical solution can
be challenging. Depending on the selected discretization and interpolation
method, numerical advection may introduce artificial diffusion, dispersive
oscillations, interpolation errors, or numerical instability.

`GeoModBox.jl` implements the following advection methods:

- upwind,
- Lax-Friedrichs (1-D only),
- staggered leapfrog,
- semi-Lagrangian, and
- passive tracer or marker-in-cell (MIC) advection.

The upwind, Lax-Friedrichs, and staggered-leapfrog methods are Eulerian
finite-difference schemes. They update scalar fields defined at the cell
centroids while the numerical grid remains fixed. In two dimensions, the
velocity components are interpolated from the staggered grid to the cell
centroids before the advection equation is solved.

The semi-Lagrangian method combines Eulerian and Lagrangian concepts. The
advected property remains defined on the fixed numerical grid, but its value
at each grid point is obtained by tracing the corresponding trajectory
backward in time to its departure point. The property is then interpolated
from the previous field at that location.

In the two-dimensional implementation of `GeoModBox.jl`, trajectories are
integrated backward using an implicit-midpoint method. The midpoint velocity
is determined iteratively from a time-centered velocity field. Velocity is
interpolated linearly during the trajectory iteration, whereas the transported
property is evaluated at the final departure points using cubic-spline
interpolation.

The passive-tracer or marker-in-cell method represents a fully Lagrangian
description of material transport. Tracers carry properties such as
temperature, phase identification, composition, or accumulated strain and
move through the fixed numerical grid according to the local velocity field.

In `GeoModBox.jl`, tracer trajectories are integrated using a classical
fourth-order Runge-Kutta method. Because the velocity components are defined
on a staggered grid, they are interpolated to the tracer positions during each
Runge-Kutta stage.

Tracer properties can subsequently be transferred back to the numerical grid.
Depending on the quantity required by the governing equations, tracer
properties may be interpolated to either cell centroids or grid vertices. For
example, density is commonly required at cell centroids, whereas viscosity may
be required at both centroids and vertices.

For further details on the implementation and usage of the different
advection methods, refer to the [1-D advection documentation](AdvOneD.md) and
the [2-D advection documentation](AdvTwoD.md).

---

## Advection Stability and Time-Step Criteria

The time-step restriction depends on the selected numerical method.

For explicit Eulerian schemes, the Courant-Friedrichs-Lewy (CFL) condition
requires that information does not travel across too many grid cells during a
single time step. In one dimension, the Courant number is defined as

$\begin{equation}
C_x = \frac{\left|v_x\right|\Delta t}{\Delta x}.
\end{equation}$

A typical stability requirement is

$\begin{equation}
C_x \leq 1.
\end{equation}$

For the two-dimensional explicit advection schemes, a practical stability
criterion is

$\begin{equation}
\frac{\max\left|v_x\right|\Delta t}{\Delta x} + \frac{\max\left|v_y\right|\Delta t}{\Delta y} \leq 1.
\end{equation}$

In the example scripts, the advection time step is commonly estimated from the
maximum velocity magnitude as

$\begin{equation}
\Delta t \leq C_{\mathrm{fac}} \frac{\min\left(\Delta x,\Delta y\right)} 
{\sqrt{ \left(\max\left|v_x\right|\right)^2 + \left(\max\left|v_y\right|\right)^2}},
\end{equation}$

where $C_{\mathrm{fac}}$ is a user-defined Courant factor.

The semi-Lagrangian method is not restricted by the conventional explicit CFL
stability condition. Nevertheless, the time step still controls the accuracy
of the calculated trajectories and the interpolation of the transported
property. Large time steps may result in inaccurate departure points or cause
them to leave the available interpolation domain.

Tracer advection using an explicit Runge-Kutta method is likewise subject to
an accuracy-related Courant restriction. Restricting the particle displacement
to approximately one grid cell per time step generally improves trajectory
accuracy and prevents tracers from crossing unresolved velocity variations.

An important consideration in all advection problems is the preservation of
the shape and amplitude of the transported field. This is particularly evident
in rigid-body rotation tests, for which the anomaly should return to its
initial position after one complete revolution without changing its shape or
amplitude. Deviations from the initial condition reveal numerical diffusion,
dispersive oscillations, trajectory errors, or interpolation errors associated
with the selected method.

---

## Examples

- [2-D advection with an analytical velocity field](../examples/Advection/Advection2D.md)
- [Resolution test of 2-D advection](../examples/Advection/AdvectionRestest2D.md)

See the [examples documentation](../examples/Examples.md) for additional
details.

---

## Exercises

- [1-D advection of a Gaussian or block anomaly](../exercises/06_1D_Advection.md)
- [2-D coupled advection-diffusion](../exercises/07_2D_Energy_Equation.md)
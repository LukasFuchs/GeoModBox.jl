# Advection Equation

In most geodynamical contexts, material is in motion, and physical properties such as temperature, density, or composition must be advected with the flow.

## Eulerian and Lagrangian Formulations

Advection describes the transport of scalar or vector quantities due to the motion of a fluid. In the following, we focus on the advection of temperature. Depending on the reference frame, the governing equation takes a different form.

### Eulerian Formulation

In a fixed (Eulerian) reference frame, the local rate of change of temperature is governed by:

$\begin{equation}
\frac{\partial T}{\partial t} = - v_i \frac{\partial{T}}{\partial{x_i}},
\end{equation}$

where $v_i$ is the velocity in the direction of $i$ and $T$ is the temperature field. Eulerian formulations are solved on a fixed grid using schemes that aim to balance accuracy, stability, and computational efficiency. `GeoModBox.jl` supports several such schemes appropriate for a range of applications.

### Lagrangian Formulation

In a moving (Lagrangian) reference frame, which follows individual fluid particles, temperature evolution is described by the material (substantial) derivative:

$\begin{equation}
\frac{DT}{Dt},
\end{equation}$

which relates to the Eulerian description via:

$\begin{equation}
\frac{DT}{Dt} = \frac{\partial T}{\partial t} + v_j \frac{\partial{T}}{\partial{x_j}}.
\end{equation}$

In the Lagrangian frame, advection reduces to solving a system of ordinary differential equations (ODEs) for particle trajectories:

$\begin{equation}
\frac{D x_j}{Dt} = v_j,
\end{equation}$

where $x_j$ is the coordinate in the direction of $j$.

## Discretization Schemes

Although simple in form, the advection equation is challenging to solve numerically. The choice of discretization and interpolation schemes can introduce numerical artifacts such as diffusion, dispersion, or instability.

For the advection, `GeoModBox.jl` includes the following numerical schemes:

- Upwind  
- Lax  
- Staggered leapfrog  
- Semi-Lagrangian  
- Passive tracer/marker method

The first four schemes work for any scalar field defined at **centroids**, including ghost nodes, and use centroid-defined velocity fields.

**Passive tracers** can be used to advect temperature or phase identifiers (IDs). When advecting phase IDs, the tracer data must also include material parameters such as viscosity or density. The physical tracer properties are interpolated either to centroids or to vertices, depending on the quantity (e.g., viscosity at vertices). The **tracer or marker-in-cell (MIC) method** tracks material properties along particle paths and solves the associated ODE system using standard time integration methods, including:

- Euler method  
- Runge–Kutta methods

The `GeoModBox.jl` focuses on the tracer advection using a fourth-order Runge–Kutta method, with velocities from the staggered grid. For more details on the tracers advection routine and how it is used in `GeoModBox.jl` refer to the [2D advection documentation](AdvTwoD.md).

An important consideration in advection is the **preservation of amplitude and shape**, especially in problems involving rigid body rotation. Numerical diffusion and interpolation artifacts can significantly affect the solution quality depending on the chosen scheme, making the selection of an appropriate method crucial.

## Advection Stability Criterion

To avoid spurious oscillations in the advected field between adjacent grid points, the *Courant–Friedrichs–Lewy (CFL)* criterion must be satisfied:

$\begin{equation}
\Delta{t} \le \frac{\textrm{min}\left(\Delta{x},\Delta{y}\right)}{\textrm{max}\left(v_x,v_y\right)},
\end{equation}$

where $\Delta{x}$ and $\Delta{y}$ are the spatial resolutions in the $x$- and $y$-directions, and $v_x$, $v_y$ are the corresponding velocity components.

This criterion ensures that no particle is advected a greater distance than the minimum grid spacing within one time step.

## Examples

- [2-D advection with constant velocity field](./examples/Advection2D.md)  
- [Resolution test of 2-D advection](./examples/AdvectionRestest2D.md)

See the [examples documentation](./Examples.md) for further details.

## Exercises

- [1-D Gaussian or block anomaly advection](../man/exercises/06_1D_Advection.md)  
- [2-D coupled advection-diffusion](../man/exercises/07_2D_Energy_Equation.md)

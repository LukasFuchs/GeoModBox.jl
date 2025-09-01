# Advection Equation

In cases where the material remains stationary (e.g., during the thermal evolution of intrusions or within a non-deforming lithosphere), it is sufficient to solve only the diffusive term of the energy equation. However, in most geodynamical contexts, material is in motion, and physical properties such as temperature, density, or composition must be advected with the flow.

## Eulerian and Lagrangian Formulations

Advection describes the transport of scalar or vector quantities due to the motion of a fluid. Depending on the reference frame, the governing equations take different forms.

### Eulerian Formulation

In a fixed (Eulerian) reference frame, the local rate of change of temperature, for example, is governed by the advection equation:

$\begin{equation}
\frac{\partial T}{\partial t} = - \vec{v} \cdot \nabla{T},
\end{equation}$

where $\vec{v}$ is the velocity vector and $T$ is the temperature field. Eulerian formulations are solved on a fixed grid using schemes that aim to balance accuracy, stability, and computational efficiency. `GeoModBox.jl` supports several such schemes appropriate for a range of applications.

### Lagrangian Formulation

In a moving (Lagrangian) reference frame, which follows individual fluid particles, temperature evolution is described by the material (substantial) derivative:

$\begin{equation}
\frac{DT}{Dt},
\end{equation}$

which relates to the Eulerian description via:

$\begin{equation}
\frac{DT}{Dt} = \frac{\partial T}{\partial t} + \vec{v} \cdot \nabla{T}.
\end{equation}$

In the Lagrangian frame, advection reduces to solving a system of ordinary differential equations (ODEs) for particle trajectories:

$\begin{equation}
\frac{Dx_i}{Dt} = v_i,
\end{equation}$

where $x_i$ are the coordinates and $v_i$ the corresponding velocity components.

## Discretization Schemes

Although simple in form, the advection equation is challenging to solve numerically. The choice of discretization and interpolation schemes—particularly when coupling grid-based fields with Lagrangian tracers—can introduce numerical artifacts such as diffusion, dispersion, or instability.

To promote clarity and modularity, `GeoModBox.jl` employs an **operator-splitting** strategy. This approach decouples the advective and diffusive terms of the temperature conservation equation and solves them sequentially. First, the advective (convective) term is solved, followed by the diffusive term. The latter is handled using the schemes described in the [Diffusion Equation documentation](./DiffMain.md). 

> **Note:** The energy equation can also be solved for diffusion and advection simultaneously using combined schemes. Interestingly, the *Forward in Time and Centered in Space (FTCS)* scheme—although unstable for pure advection—can exhibit numerical stability due to diffusion when both processes are active.

For the advection term, `GeoModBox.jl` includes the following numerical schemes:

- Upwind  
- Lax  
- Staggered leapfrog  
- Semi-Lagrangian  
- Passive tracer/marker method

The first four schemes work for any scalar field defined at **centroids**, including ghost nodes, and use centroid-defined velocity fields.

**Passive tracers** can be used to advect temperature or phase identifiers. When advecting phase IDs, the tracer data must also include material parameters such as viscosity or density. These values are interpolated either to centroids or to vertices, depending on the quantity (e.g., viscosity at vertices). The **tracer or marker-in-cell (MIC) method** tracks material properties along particle paths and solves the associated ODE system using standard time integration methods, including:

- Euler method  
- Runge–Kutta methods

The `GeoModBox.jl` focuses on the tracer advection using a fourth-order Runge–Kutta method, with velocities from the staggered grid.

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

- [1-D Gaussian or block anomaly advection](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/06_1D_Advection.ipynb)  
- [2-D coupled advection-diffusion](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/07_2D_Energy_Equation.ipynb)

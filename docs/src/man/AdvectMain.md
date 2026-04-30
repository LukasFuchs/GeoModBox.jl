# Advection Equation

In most geodynamical applications, material is in motion and physical properties such as temperature, density, or composition are transported by the flow field. This process, referred to as *advection*, describes the transport of scalar, vector, or tensor quantities due to fluid motion. In the following, we focus on the advection of temperature. However, any physical property can be advected using the methods described here. 

Depending on the chosen reference frame, the governing equation of advection takes different forms.

## Eulerian and Lagrangian Formulations

### Eulerian Formulation

In a fixed (Eulerian) reference frame, the local rate of change of temperature is governed by

$\begin{equation}
\frac{\partial T}{\partial t} = - v_i \frac{\partial{T}}{\partial{x_i}},
\end{equation}$

where $v_i$ denotes the velocity component in the $i$-direction and $T$ is the temperature field. For incompressible flow ($\nabla \cdot \mathbf{v} = 0$), the conservative and advective forms of the advection equation are equivalent. The advection equation is a first-order hyperbolic partial differential equation, and information propagates along characteristics defined by the velocity field.

Eulerian formulations are solved on a fixed numerical grid using discretization schemes designed to balance accuracy, stability, and computational efficiency. `GeoModBox.jl` provides several such schemes suitable for a wide range of applications.

### Lagrangian Formulation

In a moving (Lagrangian) reference frame, which follows individual fluid particles, temperature evolution is described by the material (substantial) derivative

$\begin{equation}
\frac{DT}{Dt},
\end{equation}$

which relates to the Eulerian description through

$\begin{equation}
\frac{DT}{Dt} = \frac{\partial T}{\partial t} + v_j \frac{\partial{T}}{\partial{x_j}}.
\end{equation}$

In the Lagrangian frame, advection reduces to solving a system of ordinary differential equations (ODEs) for the particle trajectories:

$\begin{equation}
\frac{D x_j}{Dt} = v_j,
\end{equation}$

where $x_j$ denotes the spatial coordinate in the $j$-direction.

## Discretization Schemes

Although simple in form, the advection equation is numerically challenging. The choice of discretization and interpolation schemes can introduce numerical artifacts such as artificial diffusion, dispersion, or instability.

`GeoModBox.jl` implements the following advection schemes:

- Upwind  
- Lax (1D only)
- Staggered leapfrog  
- Semi-Lagrangian  
- Passive tracer/marker-in-cell (MIC) method  

The first four schemes operate on scalar fields defined at cell centroids (including ghost nodes) and use velocities defined at the centroids. The velocities at the centroids can be interpolated from the staggered grid. 

Passive tracers can be used to advect temperature or phase identifiers (IDs). When advecting phase IDs, tracers must additionally store material properties such as viscosity or density. These tracer properties are interpolated either to cell centroids or to vertices, depending on the required quantity (e.g., viscosity evaluated at vertices).

The tracer or marker-in-cell (MIC) method tracks material properties along particle trajectories and solves the associated ODE system using standard time integration schemes, including:

- Euler method  
- Runge–Kutta methods  

In `GeoModBox.jl`, tracer advection is performed using a fourth-order Runge–Kutta scheme with velocities interpolated from the staggered grid. For further details on the tracer advection implementation and its usage, see the [2D advection documentation](AdvTwoD.md).

---

## Advection Stability Criterion

To prevent spurious oscillations between adjacent grid points, the *Courant–Friedrichs–Lewy (CFL)* condition must be satisfied:

$\begin{equation}
\Delta{t} \le \frac{\textrm{min}\left(\Delta{x},\Delta{y}\right)}{\textrm{max}\left(v_x,v_y\right)},
\end{equation}$

where $\Delta{x}$ and $\Delta{y}$ denote the grid spacing in the $x$- and $y$-directions, and $v_x$ and $v_y$ are the corresponding velocity components. This condition ensures that no particle travels farther than the minimum grid spacing within a single time step.

An important consideration in advection problems is the preservation of amplitude and shape, particularly in cases such as rigid body rotation. Numerical diffusion and interpolation errors can significantly degrade solution quality, depending on the chosen scheme. Careful selection of the advection method is therefore essential.

---

## Examples

- [2-D advection with constant velocity field](./examples/Advection2D.md)  
- [Resolution test of 2-D advection](./examples/AdvectionRestest2D.md)

See the [examples documentation](./Examples.md) for additional details.

---

## Exercises

- [1-D Gaussian or block anomaly advection](./exercises/06_1D_Advection.md)  
- [2-D coupled advection-diffusion](./exercises/07_2D_Energy_Equation.md)

# Advection Equation (2D)

In two spatial dimensions ($x$ and $y$), assuming incompressible flow, the advection equation for temperature is given by

$\begin{equation}
\frac{\partial T}{\partial t} =
-v_x\left(\frac{\partial T}{\partial x}\right)
-v_y\left(\frac{\partial T}{\partial y}\right),
\end{equation}$

where $T$ is the temperature [K], $t$ is time [s], and $v_x$ and $v_y$ are the velocity components in the $x$- and $y$-directions, respectively. In contrast to the one-dimensional case, the velocity field generally varies in both space and time.

The finite-difference discretizations of the different advection schemes closely follow their one-dimensional counterparts. For a more detailed derivation of the individual schemes and their theoretical background, the reader is referred to the [1D documentation](./AdvOneD.md).

The following sections provide an overview of the two-dimensional advection schemes currently implemented in `GeoModBox.jl`.

# Discretization Schemes

The global indexing of the central reference point $I$ follows the convention introduced in the [general solution section](./GESolution.md). The adjacent indices are defined as

$\begin{equation}
\begin{split}
I^\textrm{S} &= I^\textrm{C} - nc_x,\\
I^\textrm{W} &= I^\textrm{C} - 1,\\
I^\textrm{C} &= I,\\
I^\textrm{E} &= I^\textrm{C} + 1,\\
I^\textrm{N} &= I^\textrm{C} + nc_x,
\end{split}
\end{equation}$

where $I^\textrm{S}$, $I^\textrm{W}$, $I^\textrm{E}$, and $I^\textrm{N}$ denote the south, west, east, and north neighboring nodes of the central grid point, respectively. Together, these nodes form the five-point stencil used by the finite-difference discretizations presented below.

For the Eulerian advection schemes, the velocity at the cell centroids is obtained by bilinearly interpolating the velocity components defined on the staggered grid. These centroid velocities are subsequently used to solve the advection equation.

Tracer advection, in contrast, directly uses the staggered-grid velocity field. During the fourth-order Runge-Kutta integration, the staggered velocity components are bilinearly interpolated to the tracer positions at each Runge-Kutta stage.

## Upwind

In two dimensions, the advection equation can be discretized using the upwind scheme as

$\begin{equation}
\frac{T_{I^\textrm{C}}^{n+1}-T_{I^\textrm{C}}^n}{\Delta t} = -v_{x;I^\textrm{C}}
\begin{cases}
\frac{T_{I^\textrm{C}}^{n}-T_{I^\textrm{W}}^n}{\Delta x} &\text{if } v_{x;I^\textrm{C}} \gt 0 \\ \frac{T_{I^\textrm{E}}^{n}-T_{I^\textrm{C}}^n}{\Delta x} &\text{if } v_{x;I^\textrm{C}} \lt 0
\end{cases} 
-v_{y;I^\textrm{C}}
\begin{cases}
\frac{T_{I^\textrm{C}}^{n}-T_{I^\textrm{S}}^n}{\Delta y} &\text{if } v_{y;I^\textrm{C}} > 0 \\ 
\frac{T_{I^\textrm{N}}^{n}-T_{I^\textrm{C}}^n}{\Delta y} &\text{if } v_{y;I^\textrm{C}} < 0
\end{cases},
\end{equation}$

where $T$ denotes the temperature, $v$ the velocity, $n$ the current time step, and $\Delta t$ the time-step size. The upwind scheme is first-order accurate in both space and time and is conditionally stable. Owing to its one-sided spatial discretization, the method introduces numerical diffusion, the magnitude of which depends on both the spatial resolution and the Courant number. For implementation details, see the [source code](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/AdvectionEquation/2Dsolvers.jl).

For explicit finite-difference schemes in two dimensions, stability requires satisfaction of the Courant-Friedrichs-Lewy (CFL) condition

$\begin{equation}
|\alpha_x| + |\alpha_y| \le 1,
\end{equation}$

where

$\begin{equation}
\alpha_x = \frac{\max |v_x|\,\Delta t}{\Delta x}, \qquad 
\alpha_y = \frac{\max |v_y|\,\Delta t}{\Delta y},
\end{equation}$

and $\max |v_x|$ and $\max |v_y|$ denote the maximum absolute velocity components within the computational domain.

In the example implementations provided with `GeoModBox.jl`, the advection time step is estimated conservatively such that the maximum velocity magnitude traverses at most the minimum grid spacing during a single time step:

$\begin{equation}
\Delta t \le 
\frac{\min(\Delta x,\Delta y)}
{\sqrt{(\max |v_x|)^2 + (\max |v_y|)^2}}.
\end{equation}$

This estimate satisfies the CFL criterion for uniform Cartesian grids while providing a simple and robust time-step selection based on the maximum velocity magnitude.

## Staggered Leapfrog

The two-dimensional staggered leapfrog discretization is given by

$\begin{equation}
\frac{T_{I^\textrm{C}}^{n+1} - T_{I^\textrm{C}}^{n-1}}{2\Delta t} = 
-v_{x;I^\textrm{C}}\frac{T_{I^\textrm{E}}^{n} - T_{I^\textrm{W}}^{n}}{2\Delta x} 
-v_{y;I^\textrm{C}}\frac{T_{I^\textrm{N}}^{n} - T_{I^\textrm{S}}^{n}}{2\Delta y}.
\end{equation}$

The staggered leapfrog scheme is second-order accurate in both space and time. In contrast to the upwind scheme, it introduces very little numerical diffusion and therefore preserves the amplitude of smoothly varying fields considerably better. However, the scheme is dispersive and may produce non-physical oscillations near sharp gradients or discontinuities.

As in the one-dimensional implementation, the leapfrog scheme requires the solution at two previous time levels. Consequently, the first time step must be initialized using a separate first-order method before the leapfrog update can be applied.

For implementation details, see the [source code](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/AdvectionEquation/2Dsolvers.jl).

## Semi-Lagrangian

In two dimensions, the velocity field generally varies in both space and time. Consequently, the departure point of a characteristic trajectory cannot be determined analytically, as in the one-dimensional case. Instead, `GeoModBox.jl` employs an *implicit midpoint iteration* to determine the departure point of each characteristic trajectory. The procedure consists of the following steps:

**1. Compute a time-centered velocity field**

A temporally centered velocity field is first constructed by averaging the velocity fields at time levels $n$ and $n+1$:

$\begin{equation}
\left(v_x,v_y\right)^{n+1/2} =
\frac{\left(v_x,v_y\right)^n +
\left(v_x,v_y\right)^{n+1}}{2}.
\end{equation}$

This time-centered velocity field is used throughout the midpoint iteration.

**2. Determine the midpoint position**

Starting from the arrival point $\mathbf{x}_{I^\textrm{C}}^{n+1}$, an initial estimate of the midpoint position $X'$ of the characteristic trajectory is obtained. The velocity at $X'$ is then computed by bilinearly interpolating the time-centered velocity field to the current midpoint estimate. The midpoint position is subsequently updated using the interpolated velocity. This procedure is repeated until the midpoint position converges or until a prescribed maximum number of iterations (typically five) has been performed.

**3. Compute the departure point**

Once the midpoint iteration has converged, the departure point $\mathbf{X}^n$ at time level $n$ is computed using the converged midpoint velocity:

$\begin{equation}
\mathbf{X}^n =
\mathbf{x}_{I^\textrm{C}}^{n+1} -
\Delta t\,\mathbf{v}_{X'}^{n+1/2}.
\end{equation}$

**4. Interpolate the temperature**

Finally, the temperature field at time level $n$ is interpolated to the departure point $\mathbf{X}^n$ using cubic interpolation. The interpolated value defines the temperature at the Eulerian grid point $\mathbf{x}_{I^\textrm{C}}^{n+1}$.

Compared to explicit Eulerian schemes, the semi-Lagrangian method significantly reduces numerical diffusion and is not restricted by the conventional CFL stability criterion. However, the accuracy of the solution depends on both the trajectory integration and the interpolation of the transported field. For implementation details, see the [source code](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/AdvectionEquation/2Dsolvers.jl). 

## Passive Tracers

In two dimensions, passive tracers are advected using the classical fourth-order Runge-Kutta (RK4) method:

$\begin{equation}\begin{split}
x_p^{n+1} & = x_p^n + \frac{1}{6}k_1 + \frac{1}{3}k_2 + \frac{1}{3}k_3 + \frac{1}{6}k_4, \\
y_p^{n+1} & = y_p^n + \frac{1}{6}m_1 + \frac{1}{3}m_2 + \frac{1}{3}m_3 + \frac{1}{6}m_4,
\end{split}
\end{equation}$

where

$\begin{equation}
\begin{split}
k_1 & = \Delta t \cdot v_x(t^n,(x_p^n,y_p^n)) \\
m_1 & = \Delta t \cdot v_y(t^n,(x_p^n,y_p^n)) \\ \newline
k_2 & = \Delta t \cdot v_x(t^n+\Delta t/2,(x_p^n+k_1/2,y_p^n+m_1/2)) \\
m_2 & = \Delta t \cdot v_y(t^n+\Delta t/2,(x_p^n+k_1/2,y_p^n+m_1/2)) \\ \newline
k_3 & = \Delta t \cdot v_x(t^n+\Delta t/2,(x_p^n+k_2/2,y_p^n+m_2/2)) \\
m_3 & = \Delta t \cdot v_y(t^n+\Delta t/2,(x_p^n+k_2/2,y_p^n+m_2/2)) \\ \newline
k_4 & = \Delta t \cdot v_x(t^n+\Delta t,(x_p^n+k_3,y_p^n+m_3)) \\
m_4 & = \Delta t \cdot v_y(t^n+\Delta t,(x_p^n+k_3,y_p^n+m_3)).
\end{split}
\end{equation}$

By default, tracer trajectories are computed using the staggered-grid velocity field. The velocity components are bilinearly interpolated from the surrounding staggered grid nodes to the tracer position at each Runge-Kutta stage. Alternatively, tracer advection can also be performed using the centroid velocity field or a combination of both, depending on the application.

The tracer routines support shared-memory parallelization using Julia's multithreading capabilities. The number of threads can be specified in the VS Code Julia extension settings (`julia.NumThreads`) or by setting the corresponding Julia environment variable before starting the program.

Passive tracers can carry arbitrary material properties, such as temperature or phase identification numbers (IDs). Rather than advecting the absolute temperature, the current implementation transports the incremental temperature change between two consecutive time steps. This operator-splitting approach first computes the temperature evolution on the Eulerian grid (e.g., due to diffusion or internal heat production), after which only the temperature increment is interpolated to the tracers and advected. The updated tracer temperatures are subsequently mapped back to the computational grid. This procedure avoids repeatedly interpolating the complete temperature field and significantly reduces interpolation-induced diffusion during long-term simulations.

Phase IDs can be advected to represent compositional heterogeneity and to assign material properties such as density, viscosity, thermal conductivity, or radiogenic heat production. Depending on the governing equation being solved, tracer properties can be interpolated either to cell centroids or to grid vertices. For centroid-based quantities, the extended centroid field is used to ensure a consistent treatment of boundary cells.

When mapping tracer properties back to the finite-difference grid, several averaging schemes are available, including arithmetic, harmonic, and geometric. The appropriate averaging strategy depends on the physical property being interpolated and the governing equations in which it is used. Particular care should be taken when interpolating properties that directly influence the solution of the governing equations, such as viscosity in the momentum equation, since the chosen averaging scheme can significantly affect the numerical solution.

Currently, tracers are neither inserted nor removed dynamically during the simulation. Consequently, regions experiencing strong extension may become undersampled, whereas highly compressive regions may accumulate excessive numbers of tracers. Users should therefore ensure that a sufficiently high initial tracer density is chosen for the intended application.

For implementation details, please refer to the [source code](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/AdvectionEquation/2Dsolvers.jl).
# Advection Equation (2D)

In two spatial dimensions ($x$ and $y$), assuming incompressible flow, the advection equation for temperature is given by

$\begin{equation}
\frac{\partial{T}}{\partial{t}} = -v_x \left(\frac{\partial{T}}{\partial{x}}\right) - v_y \left(\frac{\partial{T}}{\partial{y}}\right),
\end{equation}$

where $T$ is the temperature [K], $t$ is the time [s], and $v_x$ and $v_y$ are the velocities in the $x$- and $y$-direction, respectively. In contrast to the 1D case, velocity generally varies in space and time in 2D applications.

The finite difference approximations of the different discretization schemes do not vary significantly from the 1D case. For more theoretical background, please refer to the [1D documentation](./AdvOneD.md).

The following sections provide a brief overview of the available 2D discretization schemes in the `GeoModBox.jl`. 

# Discretization Schemes

The global indexing of the central reference point $I^\textrm{C}$ used in the advection equation follows the indexing as described in the [general solution section](./GESolution.md). The indices of the adjacent points are then defined by:

$\begin{equation}\begin{split}
I^\textrm{S} & = I^\textrm{C} - nc_x,\\\
I^\textrm{W} & = I^\textrm{C} - 1,\\\   
I^\textrm{E} & = I^\textrm{C} + 1,\\\
I^\textrm{N} & = I^\textrm{C} + nc_x,
\end{split}\end{equation}$

where $I^\textrm{S}$, $I^\textrm{W}$, $I^\textrm{E}$, and $I^\textrm{N}$ are the points South, West, East, and North. These are the global indices of the five-point stencil used in the discretized FD equations below. The velocity at the centroids is computed using a bilinear interpolation of the velocities defined on the staggered grid. The staggered-grid velocity is used exclusively for tracer advection via the fourth-order Runge–Kutta scheme. For the remaining advection schemes, the centroid velocity is used. 

## Upwind

In 2D, the advection equation can be discretized using the upwind scheme as:

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

where $T$ is the temperature, $v$ the velocity, $n$ is the current time step, and $\Delta t$ the time step increment. For implementation details, see the [source code](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/AdvectionEquation/2Dsolvers.jl). The scheme is first-order accurate in space and time and is conditionally stable. It introduces numerical diffusion proportional to the grid spacing.

For explicit finite-difference schemes in two dimensions, stability requires satisfaction of the Courant–Friedrichs–Lewy (CFL) condition

$\begin{equation}
|\alpha_x| + |\alpha_y| \le 1,
\end{equation}$

with

$\begin{equation}
\alpha_x = \frac{\max |v_x|\,\Delta t}{\Delta x}, \qquad 
\alpha_y = \frac{\max |v_y|\,\Delta t}{\Delta y},
\end{equation}$

where $\max |v_x|$ and $\max |v_y|$ denote the maximum absolute velocity components over the computational domain.

In `GeoModBox.jl`, the advection time step is chosen conservatively such that the maximum velocity magnitude traverses at most the minimum grid spacing per time step:

$\begin{equation}
\Delta t \le 
\frac{\min(\Delta x,\Delta y)}
{\sqrt{(\max |v_x|)^2 + (\max |v_y|)^2}}.
\end{equation}$

## Staggered Leapfrog

The 2D staggered leapfrog discretization is given by:

$\begin{equation}
\frac{T_{I^\textrm{C}}^{n+1} - T_{I^\textrm{C}}^{n-1}}{2\Delta t} = 
-v_{x;I^\textrm{C}}\frac{T_{I^\textrm{E}}^{n} - T_{I^\textrm{W}}^{n}}{2\Delta x} 
-v_{y;I^\textrm{C}}\frac{T_{I^\textrm{N}}^{n} - T_{I^\textrm{S}}^{n}}{2\Delta y}.
\end{equation}$

For implementation details, see the [source code](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/AdvectionEquation/2Dsolvers.jl).

## Semi-Lagrangian

In 2D, velocity can vary significantly in space and time. Therefore, a modified version of the 1D semi-Lagrangian scheme is used, employing an *iterative midpoint method* to determine the origin of the characteristic trajectory. The following steps are performed:

**1. Mid-point (half time step) origin** 

Estimate the intermediate position $X'$ at time step $n+1/2$ using the velocity at the time step $n+1$ at point $I^\textrm{C}$.

**2. Mid-point velocity** 

Compute the velocity at $X'$ using temporal averaging:

$\begin{equation}
\left(v_x,v_y\right)_{I^\textrm{C}}^{n+1/2} = \frac{\left(v_x,v_y\right)_{I^\textrm{C}}^{n+1} + \left(v_x,v_y\right)_{I^\textrm{C}}^{n}}{2}.
\end{equation}$

The velocity at $X'$ can be interpolated linearly from the surrounding grid points. The procedure is repeated iteratively to update the departure point using the midpoint velocity (e.g., five iterations or until convergence). This yields an accurate temporally centered velocity.

**3. Actual origin point** 

Then, we can compute the location $X^n$ at the time step $n$ with the centered velocity: 

$\begin{equation}
\mathbf{X}^n = \mathbf{x}_{I^\textrm{C}}^{n+1} - \Delta{t}\mathbf{v}_{X'}^{n+1/2}.
\end{equation}$

**4. Update temperature** 

Use cubic interpolation to determine the temperature at the origin $X^n$, which defines the temperature at the final position $\mathbf{x}_{I^\textrm{C}}^{n+1}$.

This scheme assumes no heat sources during advection. It significantly reduces numerical diffusion but introduces interpolation-based errors. For implementation details, see the [source code](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/AdvectionEquation/2Dsolvers.jl).   

## Passive tracers

Tracers are advected in 2D using fourth-order Runge-Kutta method: 

$\begin{equation}\begin{split}
x_p^{n+1} & = x_p^n + \frac{1}{6}k_1 + \frac{1}{3}k_2 + \frac{1}{3}k_3 + \frac{1}{6}k_4, \\
y_p^{n+1} & = y_p^n + \frac{1}{6}m_1 + \frac{1}{3}m_2 + \frac{1}{3}m_3 + \frac{1}{6}m_4,
\end{split}
\end{equation}$

where:

$\begin{equation}
\begin{split}
k_1 & = \Delta{t} \cdot v_x(t^n,(x_p^n,y_p^n)) \\
m_1 & = \Delta{t} \cdot v_y(t^n,(x_p^n,y_p^n)) \\ \newline
k_2 & = \Delta{t} \cdot v_x(t^n+\Delta{t}/2,(x_p^n+k_1/2,y_p^n+m_1/2)) \\
m_2 & = \Delta{t} \cdot v_y(t^n+\Delta{t}/2,(x_p^n+k_1/2,y_p^n+m_1/2)) \\ \newline
k_3 & = \Delta{t} \cdot v_x(t^n+\Delta{t}/2,(x_p^n+k_2/2,y_p^n+m_2/2)) \\
m_3 & = \Delta{t} \cdot v_y(t^n+\Delta{t}/2,(x_p^n+k_2/2,y_p^n+m_2/2)) \\ \newline
k_4 & = \Delta{t} \cdot v_x(t^n+\Delta{t},(x_p^n+k_3,y_p^n+m_3)) \\
m_4 & = \Delta{t} \cdot v_y(t^n+\Delta{t},(x_p^n+k_3,y_p^n+m_3)) \\
\end{split}
\end{equation}$

Per default, the tracers are advected by the velocity of the staggered velocity vertices. However, if required the velocity on the centroids or a combination of both is possible. The tracer velocity is calculated using a bilinear interpolation scheme. 

The tracers can also be advected parallel by defining the maximum number of threads in the VS Code julia extension settings ("julia.NumThreads"). 

Tracers can carry different properties such as absolute temperature or phase identification number (ID). For temperature, interpolation between centroids and tracer positions is required. Note that this implementation currently does not support fully coupled temperature-momentum simulations.

Alternatively, phase IDs (linked to properties like constant density or viscosity) can be advected to simulate compositional heterogeneity. Property interpolation is performed at vertices or centroids depending on the required context (e.g., viscosity for the momentum equations). For the centroids, the extended centroid field must be used.

Caution is required when interpolating properties between the grid and tracers, especially when those properties influence the governing equations (e.g., viscosity in momentum conservation). The tracer routine enables different averaging schemes (arithmetic, harmonic, and geometric means) to map the tracer properties on the centroids or vertices using a bilinear interpolation scheme considering four surrounding finite difference grid cells for each numerical node. 

Currently, no new tracers are inserted where needed, limiting the applicability of this approach for coupled temperature-momentum models.

For implementation details, please refer to the [source code](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/AdvectionEquation/2Dsolvers.jl).
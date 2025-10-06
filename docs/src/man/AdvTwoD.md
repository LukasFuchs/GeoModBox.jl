# Advection Equation (2D)

In two dimensions ($x$ and $y$), the advection equation for the temperature conservation equation, for example, is given as follows 

$\begin{equation}
\frac{\partial{T}}{\partial{t}} = -v_x \left(\frac{\partial{T}}{\partial{x}}\right) - v_y \left(\frac{\partial{T}}{\partial{y}}\right),
\end{equation}$

where $T$ is the temperature [K], $t$ is the time [s], and $v_x$ and $v_y$ are the velocities in the $x$- and $y$-direction, respectively. 

The finite difference approximations of the different discretization schemes do not vary significantly from the 1D case. For more theoretical background, please refer to the [1D documentation](./AdvOneD.md).

The following sections provide a brief overview of the available 2D discretization schemes.

# Discretization Schemes

## Upwind

In 2D, the advection equation can be discretized using the upwind scheme as:

$\begin{equation}
\frac{T_{i,j}^{n+1}-T_{i,j}^n}{\Delta t} = -v_{x;i,j}
\begin{cases}
\frac{T_{i,j}^{n}-T_{i-1,j}^n}{\Delta x} &\text{if } v_{x;i,j} \gt 0 \\ \frac{T_{i+1,j}^{n}-T_{i,j}^n}{\Delta x} &\text{if } v_{x;i,j} \lt <0
\end{cases} 
-v_{y;i,j}
\begin{cases}
\frac{T_{i,j}^{n}-T_{i,j-1}^n}{\Delta y} &\text{if } v_{y;i,j} > 0 \\ 
\frac{T_{i,j+1}^{n}-T_{i,j}^n}{\Delta y} &\text{if } v_{y;i,j}<0
\end{cases},
\end{equation}$

where $T$ is the temperature, $v$ the velocity, $n$ is the current time step, $\Delta t$ the time step increment, and $i$, $j$ are the spatial indices in the $x$- and $y$-directions, respectively. For implementation details, see the [source code](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/AdvectionEquation/2Dsolvers.jl).

This method is stable and effective but introduces numerical diffusion if the Courant condition is not satisfied. It is also only first-order accurate in space.
   
## Staggered Leapfrog

The 2D staggered leapfrog discretization is given by:

$\begin{equation}
\frac{T_{i,j}^{n+1} - T_{i,j}^{n-1}}{2\Delta t} = 
-v_{x;i,j}\frac{T_{i+1,j}^{n} - T_{i-1,j}^{n}}{2\Delta x} 
-v_{y;i,j}\frac{T_{i,j+1}^{n} - T_{i,j-1}^{n}}{2\Delta y}.
\end{equation}$

For implementation details, see the [source code](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/AdvectionEquation/2Dsolvers.jl).

## Semi-lagragian

In 2D, velocity can vary significantly in space and time. Therefore, a modified version of the 1D semi-Lagrangian scheme is used, employing an *iterative midpoint method* to determine the origin of the characteristic trajectory. The following steps are performed:

**1. Mid-point (half time step) origin** 
Estimate the intermediate position $X'$ at time $t_{n+1/2}$ using the velocity at $t_{n+1}$ at point $(i,j)$.

**2. Mid-point velocity** 
Compute the velocity at $X'$ using temporal averaging:

$\begin{equation}
v\left(t_{n+1/2},(i,j)\right) = \frac{v\left(t_{n+1},(i,j)\right) + v\left(t_{n},(i,j)\right)}{2}.
\end{equation}$

The velocity at $X'$ can be interpolated linearly from the surrounding grid points.

**3. Actual origin point** 

Calculate the actual origin $X(t)$ for the position $x(t_{n+1},(i,j))$ using:

$\begin{equation}
X(t) = x\left(t+1,(i,j)\right) - \Delta{t}v\left(t_{n+1/2},X'\right).
\end{equation}$

This step is performed iteratively (e.g., five iterations or until convergence).

**4. Update temperature** 
Use cubic interpolation to determine the temperature at the origin $X(t)$, which defines the temperature at the final position $x(t_{n+1},(i,j))$.

This scheme assumes no heat sources during advection. It is free from numerical diffusion but introduces interpolation-based inaccuracies. For implementation details, see the [source code](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/AdvectionEquation/2Dsolvers.jl).   
## Passive tracers

Tracers are advected in 2D using fourth-oder Runge-Kutta method: 

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
m_4 & = \Delta{t} \cdot v_y(t^n+\Delta{t},(y_p^n+k_3,y_p^n+m_3)) \\
\end{split}
\end{equation}$

Per default, the tracers are advected by the velocity of the staggered velocity vertices. However, if required the velocity on the centroids or a combination of both is possible. The tracer velocity is calculated using a bilinear interpolation scheme. 

The tracers can also be advected parallel by defining the maximum number of threads in the VScode julia extension settings ("julia.NumThreads"). 

Tracers can carry different properties such as absolute temperature or phase identification number (ID). For temperature, interpolation between centroids and tracer positions is required. Note that this implementation currently does not support fully coupled temperature-momentum simulations.

Alternatively, phase IDs (linked to properties like constant density or viscosity) can be advected to simulate compositional heterogeneity. Property interpolation is performed at vertices or centroids depending on the required context (e.g., viscosity for the momentum equations).

Caution is required when interpolating properties between the grid and tracers, especially when those properties influence the governing equations (e.g., viscosity in momentum conservation).

Currently, no new tracers are inserted where needed, limiting the applicability of this approach for coupled temperature-momentum models.

For implementation details, please refer to the [source code](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/AdvectionEquation/2Dsolvers.jl).
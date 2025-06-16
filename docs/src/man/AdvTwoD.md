# Advection Equation (2D)

# Discretization Schemes

## Upwind

The idea is that the flux into the local cell will only depend on the gradient of temperature in the direction upstream. The upwind scheme is similar to a *FTCS* discretization, however, the central spacial derivatives are replaced by single-sided forward and backward finite differences and one needs to consider the advection velocity as well, to ensure that the discretization in space is always upstream. In 2-D the advection equation is then given as: 

$\begin{equation}
\frac{T_{i,j}^{n+1}-T_{i,j}^n}{\Delta t}=-v_{x;i,j}\cases{\matrix{\frac{T_{i,j}^{n}-T_{i,j-1}^n}{\Delta x} \quad \text{if} \quad v_{x;i,j} > 0 \\\ \frac{T_{i,j+1}^{n}-T_{i,j}^n}{\Delta x} \quad \text{if} \quad v_{x;i,j}<0}} 
-v_{z;i,j}\cases{\matrix{\frac{T_{i,j}^{n}-T_{i-1,j}^n}{\Delta z} \quad \text{if} \quad v_{z;i,j} > 0 \\\ \frac{T_{i+1,j}^{n}-T_{i,j}^n}{\Delta z} \quad \text{if} \quad v_{z;i,j}<0}},
\end{equation}$

where $T$ is the temperature, $v$ the velocity, $n$ is the current time step, $\Delta t$ the time step increment, and $i$ and $j$ are the indices in $z$- and $x$- direction, respectively. For more details see [*UpwindAdvection2D.m*].

This is a stable and effective way, however, with a certain amount of numerical diffusion if the *courant criteria* is not fulfilled and only first order accurate in space. The courant criteria implies that the time step is smaller than the minimum grid spacing divided by the maximum velocity, that is, a property should not be advected over a distance larger than the grid spacing, or:

$\begin{equation}
\Delta t \le \frac{\Delta x}{max(|v|)}.
\end{equation}$
   
## Staggered Leapfrog

This method considers a centered in time and centered in space discretization of the partial differentials, thus it has a higher order of accuracy in space (second order) and is suppose to not have any numerical diffusion. In 2-D the advection equation discretizes to:

$\begin{equation}
\frac{T_{i,j}^{n+1} - T_{i,j}^{n-1}}{2\Delta t}=-v_{x;i,j}\frac{T_{i,j+1}^{n} - T_{i,j-1}^{n}}{2\Delta x}-v_{z;i,j}\frac{T_{i+1,j}^{n} - T_{i-1,j}^{n}}{2\Delta z}.
\end{equation}$

For more details see [SLFAdvection2D.m].

## Semi-lagragian

This method is related to the tracer-based advection by solving ordinary differential equations (*ODEs*), where it assumes that *imaginary tracers* are located at certain positions and land directly at the finite difference grid nodes after advection within one time step. Thus, one needs to calculate the *origin points* for each grid node back in time (e.g., one Euler time step) with a given velocity field (e.g., using an *iterative mid-point scheme*, i.e. one uses the velocity at a point half a time step backward in time) and then to interpolate the property from the regular grid points to the determined *origin points*. This scheme assumes that no heat-sources were active during the advection. The method does not have any numerical diffusion but shows inaccuracies due to the interpolation method. For more details see [SemiLagAdvection2D.m].

- central point iteration
   
## Passive tracers

Here, one assumes that the model domain is completely filled with so-called *tracers* or *markers*. These tracers are then advected by solving the *ODE* of a particle advection using a certain method (e.g., Euler or Runge Kutta) and they transport any property stored on them. However, care needs to be taken when interpolating those properties from the regular grid onto the tracers and back. This is even more complex if the property advected does have an effect on parameters controlling the governing equations (e.g., the viscosity in continuum euqation).

Here, I advect the tracers using Runge-Kutta fourth order; the tracers can transport the absolute temperature and the composition (so far only for two compositions with a constant viscosity and density). The property is then interpolated back to the regular grid points every time step. For more details see [AdvectMarker2D.m] and [TracerInterp.m].
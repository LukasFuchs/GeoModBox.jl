# Energy Conservation Equation

The conservation of energy is a fundamental physical principle stating that energy can neither be created nor destroyed, but only transformed. In geodynamical modeling, this principle is typically expressed in terms of **temperature**, which evolves due to **diffusive** and **advective** transport processes. Assuming radioactive heat sources only, the general energy equation is defined as:

$\begin{equation}
\left(\frac{\partial E}{\partial t} + v_j\frac{\partial{E}}{\partial{x_j}}\right) + \frac{\partial q_{i}}{\partial x_{i}} = \rho H,
\end{equation}$

where energy is defined as $E = c_p \rho T$. Here 
$c_p$ is the specific heat capacity [J/kg/K],
$\rho$ is the density [kg/m³],
$T$ is temperature [K],
$\partial/\partial t$ is the time derivative,
$t$ is time [s],
$v_j$ is the velocity [m/s] in direction $j$,
$q_i$ is the heat flux [W/m²] in direction $i$,
$\partial/\partial x_i$ is a directional derivative in $i$, and
$H$ is the internal heat production per unit mass [W/kg]. Repeated indices imply summation.

The heat flux $q_i$ is described by **Fourier’s law**:

$\begin{equation}
q_i = - k \frac{\partial{T}}{\partial{x_i}},
\end{equation}$

where $k$ is the thermal conductivity [W/m/K]. The flux is directed opposite to the temperature gradient and represents the amount of heat passing through a unit surface per unit time.

Substituting Fourier’s law into the energy equation yields the **temperature conservation equation** in Eulerian form:

$\begin{equation}
\rho c_p \left(\frac{\partial T}{\partial t} + v_j\frac{\partial{T}}{\partial{x_j}}\right) = \frac{\partial}{\partial x_i}\left(k\frac{\partial{T}}{\partial{x_i}}\right) + \rho H.
\end{equation}$

This equation captures temperature changes due to **diffusion** (right-hand side) and **advection** (left-hand side). For simplicity and assuming a spatially constant internal heat production, these processes can be split using an *operator splitting* technique, solving the advection and diffusion steps sequentially. If internal heat production varies spatially, a more advanced advection scheme is required to account for source term integration.

# Heat Diffusion Equation

In cases where the material remains stationary (e.g., during the thermal evolution of intrusions or within a non-deforming lithosphere), it is sufficient to solve only the diffusive term of the energy equation. Neglecting the advection part of the temperature conservation equation, the heat diffusion equation is defined as: 

$\begin{equation}
\rho c_p \frac{\partial T}{\partial t} = \frac{\partial}{\partial x_i}\left(k\frac{\partial{T}}{\partial{x_i}}\right) + \rho H.
\end{equation}$

```GeoModBox.jl``` provides several finite difference (FD) schemes to solve the heat diffusion equation of the time-dependent or steady-state temperature conservation equation—including optional radioactive heating and variable thermal parameters—in both [1-D](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/HeatEquation/1Dsolvers.jl) and [2-D](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/HeatEquation/2Dsolvers.jl). Available methods include:

- Forward Euler  
- Backward Euler  
- Crank–Nicolson  
- Alternating Direction Implicit (ADI)

Except for ADI, all solvers are implemented for constant and variable thermal parameters to solve a linear problem using a left-matrix division and a non-linear problem using the defection correction. See the documentation for the [1-D](./DiffOneD.md) and [2-D](./DiffTwoD.md) solvers for detailed descriptions of each method. Currently, only *Dirichlet* and *Neumann* boundary conditions are supported. See the [HeatEquation source directory](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/HeatEquation/) for implementation details.

## Examples

- [1-D oceanic geotherm](./examples/OceanicGeotherm.md)  
- [1-D continental geotherm](./examples/ContinentalGeotherm.md)  
- [Comparison of FD schemes on a Gaussian anomaly](./examples/GaussianDiffusion1D.md)  
- [2-D resolution test with Gaussian anomaly](./examples/GaussianDiffusion2D.md)  
- [2-D Poisson equation resolution test](./examples/PoissonRestest.md)

For more details, see the full [example directory](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/DiffusionEquation/).

## Exercises

- [1-D forward Euler](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/02_1D_Heat_explicit.ipynb)  
- [1-D backward Euler](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/03_1D_Heat_implicit.ipynb)  
- [2-D Poisson equation](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/04_2D_Diffusion_Stationary.ipynb)  
- [2-D transient plume heating](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/05_2D_Diffusion_TD_Plume.ipynb)  
- [2-D transient sill heating](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/05_2D_Diffusion_TD_Sill.ipynb)

# Heat Advection Equation

See the [advection documentation](./AdvectMain.md) for more theoretical and associated source codes for method-specific details. To solve the **advective** component of the temperature conservation equation, ```GeoModBox.jl``` offers several schemes:

- Upwind scheme  
- Staggered leapfrog scheme  
- Semi-Lagrangian scheme  
- Passive tracers/markers

Tracer-related functionality is located in [src/Tracers](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/Tracers/), whereas other advection schemes are implemented in [src/AdvectionEquation](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/AdvectionEquation/).

## Examples

- [2-D advection with constant velocity field](./examples/Advection2D.md)  
- [Resolution test of 2-D advection](./examples/AdvectionRestest2D.md)


## Exercises

- [1-D Gaussian or block anomaly advection](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/06_1D_Advection.ipynb)  
- [2-D coupled advection-diffusion](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/07_2D_Energy_Equation.ipynb)

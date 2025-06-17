# Energy Conservation Equation

The conservation of energy is a fundamental physical principle, stating that energy cannot be created or destroyed, only transformed. In geodynamical modeling, this is commonly expressed in terms of **temperature**, which is transported through **diffusive** and **convective** processes. Assuming radioactive heat sources only, the general energy equation is defined as:

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
$\partial/\partial x_i$ is a directional derivative in $i$,
$H$ is the internal heat production per unit mass [W/kg]. Repeated indices imply summation.

The heat flux $q_i$ is described by **Fourier’s law**:

$\begin{equation}
q_i = - k \frac{\partial{T}}{\partial{x_i}},
\end{equation}$

where $k$ is the thermal conductivity [W/m/K]. The flux is directed opposite to the temperature gradient and represents the amount of heat passing through a unit surface per unit time.

Substituting Fourier’s law into the energy equation yields the **temperature conservation equation** in Eulerian form:

$\begin{equation}
\rho c_p \left(\frac{\partial T}{\partial t} + v_j\frac{\partial{T}}{\partial{x_j}}\right) = -\frac{\partial q_i}{\partial x_i} + \rho H.
\end{equation}$

This equation captures temperature changes due to **diffusion** (right-hand side) and **advection** (left-hand side). For simplicity and assuming a spatially constant internal heat production, these processes can be split using an *operator splitting* technique, solving the advection and diffusion steps sequentially. If the internally heat production does vary one needs to consider a slightly improved advection scheme by integrating the source terms.

## Heat Diffusion Equation

```GeoModBox.jl``` provides several finite difference (FD) schemes to solve the diffusive component of the time-dependent or steady-state temperature equation (with optional radioactive heating) in both [1-D](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/HeatEquation/1Dsolvers.jl) and [2-D](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/HeatEquation/2Dsolvers.jl). Available methods include:

- Forward Euler  
- Backward Euler  
- Crank–Nicolson  
- Alternating Direction Implicit (ADI)

See the documentation for the [1-D](./DiffOneD.md) and [2-D](./DiffTwoD.md) solvers for detailed descriptions of each method.

At present, only *Dirichlet* and *Neumann* boundary conditions are supported. Most implementations assume constant thermal properties, with exceptions in some 1-D and 2-D solvers. See the [HeatEquation source directory](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/HeatEquation/) for implementation details.

### Examples

- [1-D oceanic geotherm](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/DiffusionEquation/1D/OceanicGeotherm_1D.jl)  
- [1-D continental geotherm](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/DiffusionEquation/1D/ContinentalGeotherm_1D.jl)  
- [Comparison of FD schemes on a Gaussian anomaly](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/DiffusionEquation/1D/Heat_1D_discretization.jl)  
- [2-D resolution test with Gaussian anomaly](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/DiffusionEquation/2D/Gaussian_Diffusion.jl)  
- [2-D Poisson equation resolution test](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/DiffusionEquation/2D/Poisson_ResTest.jl)

For explanations, see the [examples documentation](./Examples.md) and the full [example directory](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/DiffusionEquation/).

### Exercises

- [1-D forward Euler](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/02_1D_Heat_explicit.ipynb)  
- [1-D backward Euler](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/03_1D_Heat_implicit.ipynb)  
- [2-D Poisson equation](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/04_2D_Diffusion_Stationary.ipynb)  
- [2-D transient plume heating](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/05_2D_Diffusion_TD_Plume.ipynb)  
- [2-D transient sill heating](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/05_2D_Diffusion_TD_Sill.ipynb)

## Heat Advection Equation

To solve the **advective** component of the temperature equation, ```GeoModBox.jl``` offers several schemes:

- Upwind scheme  
- Staggered leapfrog scheme  
- Semi-Lagrangian scheme  
- Passive tracers/markers

See the [advection documentation](./AdvectMain.md) and source code for method-specific details. Tracer code resides in [src/Tracers](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/Tracers/), while the other schemes are in [src/AdvectionEquation](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/AdvectionEquation/).

### Examples

- [2-D advection with constant velocity field](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/AdvectionEquation/2D_Advection.jl)  
- [Resolution test of 2-D advection](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/AdvectionEquation/2D_Advection_ResolutionTest.jl)

See the [examples documentation](./Examples.md) for further details.

### Exercises

- [1-D Gaussian or block anomaly advection](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/06_1D_Advection.ipynb)  
- [2-D coupled advection-diffusion](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/07_2D_Energy_Equation.ipynb)

# Energy Conservation Equation

The conservation of energy is a fundamental physical principle stating that energy can neither be created nor destroyed, but only transformed. In geodynamical modeling, this principle is typically expressed in terms of **temperature**, which evolves due to **diffusive** and **advective** transport processes.

Assuming radioactive heat production as the only internal heat source, the general energy equation can be written as

$\begin{equation}
\left(\frac{\partial E}{\partial t} + v_j\frac{\partial{E}}{\partial{x_j}}\right) + \frac{\partial q_{i}}{\partial x_{i}} = \rho H,
\end{equation}$

where the internal energy is defined as $E = c_p \rho T$. Here

* $c_p$ is the specific heat capacity [J/kg/K],
* $\rho$ is the density [kg/m³],
* $T$ is the temperature [K],
* $t$ is time [s],
* $v_j$ is the velocity [m/s] in direction $j$,
* $q_i$ is the heat flux [W/m²] in direction $i$,
* $\partial/\partial x_i$ denotes the spatial derivative in direction $i$, and
* $H$ is the internal heat production per unit mass [W/kg].

Indices $i$ and $j$ denote spatial directions. Repeated indices follow Einstein summation notation, implying summation over spatial dimensions.

The heat flux $q_i$ is described by Fourier’s law

$\begin{equation}
q_i = - k_i \frac{\partial{T}}{\partial{x_i}},
\end{equation}$

where $k_i$ is the thermal conductivity [W/m/K] in direction $i$. The heat flux is directed opposite to the temperature gradient and represents the amount of heat passing through a unit surface per unit time.

Substituting Fourier’s law into the energy equation yields the **temperature conservation equation** in Eulerian form:

$\begin{equation}
\rho c_p \left(\frac{\partial T}{\partial t} + v_j\frac{\partial{T}}{\partial{x_j}}\right) = \frac{\partial}{\partial x_i}\left(k_i\frac{\partial{T}}{\partial{x_i}}\right) + \rho H.
\end{equation}$

This equation describes temperature changes due to diffusion (right-hand side) and advection (left-hand side).

In `GeoModBox.jl`, these equations are solved using finite difference discretizations on structured grids. Temperature is defined at cell centers (centroids), while fluxes are evaluated at cell interfaces using a staggered-grid formulation. This conservative discretization ensures consistent heat flux calculations and numerical stability.

For numerical implementation it is often convenient to split these processes using an operator splitting approach, in which the advection and diffusion steps are solved sequentially. Operator splitting is commonly used when diffusion and advection operate on different time scales. Spatially variable heat production may require careful treatment within the advection step. In practice, `GeoModBox.jl` advances the temperature field in time by sequentially solving the advection and diffusion components of the temperature equation using operator splitting.

> **Note:** The energy equation can also be solved for diffusion and advection simultaneously using combined schemes. Interestingly, the *Forward in Time and Centered in Space (FTCS)* scheme, although unstable for pure advection, may become numerically stable when diffusion is included.

# Heat Diffusion Equation

In situations where the material is stationary (for example during the thermal evolution of magmatic intrusions or within a non-deforming lithosphere), only the diffusive component of the energy equation needs to be solved.

Neglecting the advective term, the heat diffusion equation becomes

$\begin{equation}
\rho c_p \frac{\partial T}{\partial t} = \frac{\partial}{\partial x_i}\left(k_i\frac{\partial{T}}{\partial{x_i}}\right) + \rho H.
\end{equation}$

`GeoModBox.jl` provides several finite difference (FD) schemes to solve the heat diffusion equation for both time-dependent and steady-state problems. The implementation supports optional radioactive heating and variable thermal parameters in both

- [1-D](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/HeatEquation/1Dsolvers.jl)
- [2-D](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/HeatEquation/2Dsolvers.jl)

Available discretization methods include

- Forward Euler  
- Backward Euler  
- Crank–Nicolson  
- Alternating Direction Implicit (ADI)

Except for ADI, all solvers are implemented for both constant and variable thermal properties. Linear problems can be solved directly using a left-matrix division, while non-linear problems are solved iteratively using the defect correction method.

Detailed descriptions of the numerical schemes are provided in the documentation of the

- [1-D solvers](./DiffOneD.md)
- [2-D solvers](./DiffTwoD.md)

Currently, Dirichlet and Neumann boundary conditions are supported. Implementation details can be found in the [HeatEquation source directory](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/HeatEquation/).

## Examples

The following example scripts demonstrate the application of the diffusion solvers:

- [1-D oceanic geotherm](./examples/OceanicGeotherm.md)
- [1-D continental geotherm](./examples/ContinentalGeotherm.md)
- [Comparison of FD schemes on a Gaussian anomaly](./examples/GaussianDiffusion1D.md)
- [2-D resolution test with Gaussian anomaly](./examples/GaussianDiffusion2D.md)
- [2-D Poisson equation resolution test](./examples/PoissonRestest.md)

Additional examples can be found in the full  
[example directory](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/DiffusionEquation/).

## Exercises

The following exercises provide hands-on practice with the implemented solvers:

- [1-D forward Euler](./exercises/02_1D_Heat_explicit.md)
- [1-D backward Euler](./exercises/03_1D_Heat_implicit.md)
- [2-D Poisson equation](./exercises/04_2D_Diffusion_Stationary.md)
- [2-D transient plume heating](./exercises/05_2D_Diffusion_TD_Plume.md)
- [2-D transient sill heating](./exercises/05_2D_Diffusion_TD_Sill.md)

# Heat Advection Equation

For the **advective** component of the temperature conservation equation, `GeoModBox.jl` provides several numerical schemes. For additional theoretical background and implementation details, see the [advection documentation](./AdvectMain.md).

Available methods include

- Upwind scheme
- Staggered leapfrog scheme
- Semi-Lagrangian scheme
- Passive tracers/markers

Tracer-related functionality is implemented in

- [src/Tracers](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/Tracers/)

Other advection schemes are implemented in

- [src/AdvectionEquation](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/AdvectionEquation/)

---

## Examples

- [2-D advection with constant velocity field](./examples/Advection2D.md)
- [Resolution test of 2-D advection](./examples/AdvectionRestest2D.md)


## Exercises

- [1-D Gaussian or block anomaly advection](./exercises/06_1D_Advection.md)
- [2-D coupled advection–diffusion](./exercises/07_2D_Energy_Equation.md)
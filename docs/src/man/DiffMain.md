# Energy Conservation Equation

The conservation of energy is a fundamental principle in physics and defines that the loss and generation of energy needs to be equal within a closed system. In terms of a geodynamical problem, energy is described by temperature, which is transported mainly through *conductive* and *convective* processes, such that a general energy equation is defined as followed (assuming only radioactive heat sources):

$\begin{equation}
\left(\frac{\partial E}{\partial t} + \overrightharpoon{v} \cdot \overrightharpoon{\nabla} E\right) + \frac{\partial q_{i}}{\partial x_{i}} = \rho H,
\end{equation}$

where the energy is described as $E=c_{p} \rho T$, and $c_{p}$ is the specific heat capacity [ $J/kg/K$ ], $\rho$ is a reference density [ $kg/m^{3}$ ], $T$ is the temperature [ $K$ ], $t$ is the time [ $s$ ], $\overrightharpoon{v}$ is the velocity vector [ $m/s$ ], $q_{i}$ is the heat flux in direction of $i$ [ $W/m^{2}$ ], $\frac{\partial}{\partial{x_i}}$ is a directional derivative in direction of $i$, and $H$ the heat production rate per mass [ $W/kg$ ]. The repeated index means a summation of derivatives. This conservation law contains the variation of the heat flux in a certain direction, where the heat flux is defined by the Fourier’s law as followed: 

$\begin{equation}
\overrightharpoon{q} = - k \overrightharpoon{\nabla} T,
\end{equation}$

where $k$ is the thermal conductivity [ $W/m/K$ ]. The heat flux is the amount of heat that passes through a unit surface area, per unit time and is positive in the direction of decreasing temperature, that is in the case when the temperature gradient is negative. The *temperature conservation equation*, or *temperature equation*, in an Eulerian form can then be written as: 

$\begin{equation}
\rho c_p \left(\frac{\partial T}{\partial t} + \overrightharpoon{v} \cdot \overrightharpoon{\nabla} T\right) = -\frac{\partial q_i}{\partial x_i} + \rho H.
\end{equation}$

This form of the *temperature equation* describes the variation of temperature due to a *conductive* (right hand side of the equation) and *convective* (left hand side of the equation) process. For a matter of simplicity, one can consider those terms in a separate manner and solve the *temperature equation* using an *operator splitting* method, that is one first solves the *convective* part, followed by the *conductive* part.

## Heat Diffusion Equation

The ```GeoModBox.jl``` provides different finite difference (**FD**) schemes to solve the *diffusive part* of the time-dependent or steady-state *temperature equation* including radioactive heating, in [1-D](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/src/HeatEquation/1Dsolvers.jl) and [2-D](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/src/HeatEquation/2Dsolvers.jl). The solvers include: 

- forward Euler, 
- backward Euler, 
- Crank-Nicholson approach, and 
- Alternate Direct Implicit (ADI).

For more details regarding each FD-scheme see the documentation of the [1-D](./DiffOneD.md) and [2-D](./DiffTwoD.md) solvers, respectively. Currently, only *Dirichlet* and *Neumann* thermal boundary conditions are available and most functions assume constant thermal parameters (except for the 1-D and some 2-D solvers). For more details on the implementation of the solvers see the [src/HeatEquation](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/src/HeatEquation/) directory. 

The examples of solving the *heat diffusion equation* include, amongst others: 
- the determination of an [oceanic](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/examples/DiffusionEquation/1D/OceanicGeotherm_1D.jl) and [continental](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/examples/DiffusionEquation/1D/ContinentalGeotherm_1D.jl) 1-D geotherm profile, 
- [a comparison of the different **FD**-schemes applied on a 1-D gaussian temperature anomaly](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/examples/DiffusionEquation/1D/Heat_1D_discretization.jl), 
- [a 2-D resolution test for each **FD**-scheme using a gaussian temperature anomaly](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/examples/DiffusionEquation/2D/Gaussian_Diffusion.jl), and
- [a resolution test for a 2-D poisson problem](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/examples/DiffusionEquation/2D/Poisson_ResTest.jl). 

For detailed information and explenations of the examples see the [documentation](./Examples.md) and for additional examples see the [example directory](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/examples/DiffusionEquation/).

The [exercises](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/exercises/) include
- [the forward](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/exercises/02_1D_Heat_explicit.ipynb) and [backward](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/exercises/03_1D_Heat_implicit.ipynb) Euler method to solve the 1-D diffusion equation, 
- [a 2-D poisson problem](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/exercises/04_2D_Diffusion_Stationary.ipynb), and
- a time-dependent temperature distribution within the lithosphere assuming a [plume](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/exercises/05_2D_Diffusion_TD_Plume.ipynb) or [sill](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/exercises/05_2D_Diffusion_TD_Sill.ipynb).

## Heat Advection Equation

To solve the *advective part* of the *temperature equation*, the ```GeoModBox.jl``` provides the following different methods: 
- the upwind scheme,
- the staggered-leaped frog scheme, 
- the semi-lagrangian advection scheme, and
- passive tracers/markers. 

For more details regarding each advection scheme see the [documentation](./AdvectMain.md). The routines are structured in such a way, that any property can be advected with the first **three** advection methods listed above, as long as the property is defined on the *centroids* including *ghost nodes* on all boundaries. The velocity to advect the property also needs to be defined at the *centroids*. Using passive tracers, one can, so far, choose to either advect the (absolut) temperature or the phase ID. In case of advecting the phase ID, one must define a certain rheology ($\eta$) and/or density ($\rho$) associated to each phase. The phase ID is used to interpolate the corresponding property from the tracers to the *centroids* or the *vertices* (e.g., in the case of the viscosity). The tracers can be initialized with a certain random perturbation for their initial position and are advected using Runge-Kutta 4th order. For the tracer advection one can use the staggered velocity grid. 

For more details on the implementation of the tracer advection see [here](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/src/Tracers/2Dsolvers.jl) or [here](./AdvectMain.md). The solvers for the tracer advection method are located in [src/Tracers](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/src/Tracers/), where the remaining advection routines are located in [src/AdvectionEquation](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/src/AdvectionEquation/). 

A key aspect for the advection is the **amplitude** and the **shape** of an advected property. For example, in case of a rigid body rotation, none of them should change. Depending on the method, however, numerical diffusion or interpolation effects can lead to strong deviations. Thus, care needs to be taken with respect to solving the *advection equation*. 

The ```GeoModBox.jl``` contains several [routines](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/src/InitialCondition/2Dini.jl) to setup a certain initial anomaly, either for properties defined on their correspondig grid (i.e., temperature, velocity, or phase) or for tracers. Within the examples and the exercise one can choose different initial temperature and velocity conditions. See [here](./Ini.md) for a more detailed explenation of the initial condition setup. 

The examples for a two dimensional advection problem include:
- [a 2-D advection, assuming a constant velocity field (e.g., a rigid body rotation)](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/examples/AdvectionEquation/2D_Advection.jl), and
- [a resolution test of the same advection example](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/examples/AdvectionEquation/2D_Advection_ResolutionTest.jl). 

For detailed information and explenations of the examples see the [documentation](./Examples.md).

The exercises include a
- [1-D advection of a gaussian or block anomaly](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/exercises/06_1D_Advection.ipynb), and
- [a 2-D advection coupled with the solution of the diffusion equation](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/exercises/07_2D_Energy_Equation.ipynb).
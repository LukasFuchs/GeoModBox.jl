# GeoModBox.jl

The **Geod**ynamic **Mod**elling Tool**Box** is a julia package mainly used for teaching purposes. The package provides different finite difference, staggered, discretization schemes to numerically solve the governing equations for a two-dimensional geodynamic problem. The governing equations are the conservation equations of 
1) [**energy**](./man/DiffMain.md), 
2) [**momentum**](./man/MomentumMain.md), 
3)  [**mass** and **compositon**](./man/AdvectMain.md). 

The ```GeoModBox.jl``` includes a series of [exercises](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/exercises/) and [examples](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/examples/) of different geodynamically well defined problems. The exercises are given as Jupyter notebooks for the students to complete. The theoretical background is mainly given here in the documentation.

------------------

## Staggered Finite Difference 

------------------

## [Energy Conservation Equation](./man/DiffMain.md)

In geodynamics, the energy is described by the temperature and needs to be conserved within a closed system. Here, we solve the *temperature conservation equation*, or *temperature equation*, using an operator splitting method, that is, we first solve the *advective* part of the *temperature equation*, followed by the *diffusive* part. 

------------------

### [Heat Diffusion Equation](./man/DiffOneD.md)

The ```GeoModBox.jl``` provides different finite difference (**FD**) schemes (e.g., [forward](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/src/HeatEquation/ForwardEuler.jl) and [backward](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/src/HeatEquation/BackwardEuler.jl) Euler, [Crank-Nicholson approach](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/src/HeatEquation/CNA.jl), [ADI](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/src/HeatEquation/ADI.jl)) to solve the *diffusive part* of the time-dependent or steady-state *temperature equation* including radioactive heating, in [1-D](./man/DiffOneD.md) and [2-D](./man/DiffTwoD.md). The solvers are located in [src/HeatEquation](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/src/HeatEquation/). So far, only *Dirichlet* and *Neumann* thermal boundary conditions are available. Currently, most of the functions assume constant thermal parameters (except for the 1-D solvers). 

The examples of solving the *heat diffusion equation* include, amongst others: 
- the determination of an [oceanic](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/examples/DiffusionEquation/1D/OceanicGeotherm_1D.jl) and [continental](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/examples/DiffusionEquation/1D/ContinentalGeotherm_1D.jl) 1-D geotherm profile, 
- [a comparison of the different **FD**-schemes applied on a 1-D gaussian temperature anomaly](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/examples/DiffusionEquation/1D/Heat_1D_discretization.jl), 
- [a 2-D resolution test for each **FD**-scheme using a gaussian temperature anomaly](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/examples/DiffusionEquation/2D/Gaussian_Diffusion.jl), and
- [a resolution test for a 2-D poisson problem](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/examples/DiffusionEquation/2D/Poisson_ResTest.jl). 

For more examples see the [example folder](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/examples/DiffusionEquation/). 

The [exercises](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/exercises/) include solving 
- the 1-D diffusion equation using the [forward](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/exercises/02_1D_Heat_explicit.ipynb) and [backward](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/exercises/03_1D_Heat_implicit.ipynb) Euler methods, 
- [a 2-D poisson problem](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/exercises/04_2D_Diffusion_Stationary.ipynb), and
- a time-dependent temperature distribution within the lithosphere assuming a [plume](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/exercises/05_2D_Diffusion_TD_Plume.ipynb) or [sill](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/exercises/05_2D_Diffusion_TD_Sill.ipynb).

### [Heat Advection Equation](./man/AdvectMain.md)

To solve the *advective part* of the *temperature equation*, the ```GeoModBox.jl``` provides the following different methods: 
- an upwind scheme,
- the staggered -leaped frog scheme, 
- a semi-lagrangian advection, and
- passive tracers/markers. 

The solvers for a the tracer advection method are located in [src/Tracers](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/src/Tracers/), where the remaining advection routines are located in [src/AdvectionEquation](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/src/AdvectionEquation/). The routines are structured in such a way, that any property, as long as the property is defined on the *centroids* including *ghost nodes* on all boundaries, can be advected with the first **three** advection methods listed above. Using passive tracers, one can, so far, choose to either advect the temperature or phases. In case of advecting phases, one can define a certain rheology ($\eta$) or density ($\rho$) associated to each phase. The phase ID is used to interpolate the corresponding property from the tracers to the *centroids*. 

A key aspect for the advection equation is the conservation of the **amplitude** and the **shape** of an anomaly. Depending on the method, numerical diffusion or interpolation effects lead to strong deviations of the initial anomaly. For more details see [here](./man/AdvectMain.md). The ```GeoModBox.jl``` contains several [routines](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/src/InitialCondition/2Dini.jl) to setup a certain initial anomaly, either for properties defined on their correspondig grid (i.e., temperature, velocity, or phase) or for tracers advecting (so far) a certain temperature or phase. Within the examples and the exercise one can choose different initial temperature and velocity conditions.

The examples for a two dimensional advection problem include:
- [a 2-D advection, assuming a constant velocity field (e.g., a rigid body rotation)](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/examples/AdvectionEquation/2D_Advection.jl), and
- [a resolution test of the same advection example](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/examples/AdvectionEquation/2D_Advection_ResolutionTest.jl). 

The exercises include a: 
- [1-D advection of a gaussian or block anomaly](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/exercises/06_1D_Advection.ipynb), and
- [a 2-D advection coupled with the solution of the diffusion equation](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/exercises/07_2D_Energy_Equation.ipynb).

------------------

## [Momentum Conservation Equation](./man/MomentumMain.md)

------------------

## Code Structure

### Initial Conditions 

### Scaling

------------------
------------------
## [Benchmarks and Examples](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/examples/)

### [Gaussian Temperature Diffusion](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/examples/DiffusionEquation/2D/Gaussian_Diffusion.jl)

![GaussianDiffusion](./assets/Gaussian_Diffusion_CNA_nx_100_ny_100.gif)

**Figure 1. Gaussian Diffusion.** Time-dependent, diffusive solution of a 2-D Gaussian temperature anomaly using the [Crank-Nicholson approach](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/src/HeatEquation/CNA.jl) in comparison to its analytical solution. Top Left: 2-D temperature field of the numerical solution and isotherms lines of the numerical (solid black) and analytical (dashed yellow) solution. Top Right: Total deviation to the analytical solution. Bottom Left: 1-D y-profile along x=0. Bottom Right: Root Mean Square total devation of the temperature over time. 

![GDResTest](./assets/Gaussian_ResTest.png)

**Figure 2. Resolution test.** Maximum *RMS* $\varepsilon$, maximum, and mean temperature for each **FD**-scheme and multiple resolutions. 

### [Rigid-Body-Rotation](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/examples/AdvectionEquation/2D_Advection.jl)

![RigidBodyI](./assets/2D_advection_circle_RigidBody_upwind_100_100_nth_1.gif)

![RigidBodyII](./assets/2D_advection_circle_RigidBody_semilag_100_100_nth_1.gif)

![RigidBodyIII](./assets/2D_advection_circle_RigidBody_tracers_100_100_nth_1.gif)

**Figure 3. Rigid-Body-Rotation.** Time-dependent solution of a rotating circular temperature anomaly using the **upwind (first)**, **semi-lagrangian (second)**, and **tracer (third)** method. Within a circular area of our model domain the velocity is set to the velocity of a rigid rotation and outside euqal to zero. The temperature is scaled by the maximum temperature of the anomaly. 

### [Falling Block](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/examples/StokesEquation/2D/FallingBlockBenchmark.jl)

![FallingBlockTD](./assets/Falling_block_ηr_0.0_tracers.gif)

**Figure 4. Isoviscous Falling Block.** Time-dependent solution of an isoviscous falling block example. The problem is solved with a solver for variable viscosities. The tracers advect the phase ID, which is used to interpolate the density and viscosity on the centroids and vertices, respectively. 

![FBSinkinVeloc](./assets/FallingBlock_SinkingVeloc_tracers.png)

**Figure 5. Falling Block Sinking Velocity.** Sinking velocity of the block with respect to the viscosity ratio $\eta_r$ at the initial condition. 

![FBFinalStage](./assets/FallingBlock_FinalStage_tracers.png)

**Figure 6. Falling Block Benchmark.** Final tracers distribution for specific cases with $\eta_r \ge 0$. 

------------------
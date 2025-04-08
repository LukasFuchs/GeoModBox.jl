# GeoModBox.jl

The **Geod**ynamic **Mod**elling Tool**Box** is a julia package mainly used for teaching purposes. The package provides different finite difference, staggered, discretization schemes to numerically solve the governing equations for a two-dimensional geodynamic problem. The governing equations are the conservation equations of 

1) [**energy**](./man/DiffMain.md), 
2) [**momentum**](./man/MomentumMain.md), 
3)  [**mass** and **compositon**](./man/AdvectMain.md). 

The ```GeoModBox.jl``` includes a series of [exercises](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/) and [examples](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/) of different geodynamically well defined problems. The exercises are given as Jupyter notebooks for the students to complete. The theoretical background is mainly given here in the documentation.

## [Staggered Finite Difference](./man/GESolution.md)

To properly solve the governing equations, a staggered finite difference scheme is choosen for the *energy* and *momentum* equations. A staggered grid enables a correct, straight forward implementation of certain boundary conditions and enables the conservation of stress between the nodes in case of a variable viscosity. This also requires that certain parameters are defined on different grids. 

Here, the temperature, density, pressure are defined on the *centroids*, ... are defined on the *vertices*, the velocities are defined in between the *vertices*, and the viscosity is needed on both. 

For more details on how this is used in the ```GeoModBox.jl``` see [here](./man/GESolution.md).

## Energy Conservation Equation

In geodynamics, the energy is described by the temperature and needs to be conserved within a closed system. Here, we solve the *temperature conservation equation*, or *temperature equation*, using an *operator splitting* method, that is, we first solve the *advective* part of the *temperature equation*, followed by the *diffusive* part. 

### [Heat Diffusion Equation](./man/DiffMain.md)

The ```GeoModBox.jl``` provides different finite difference (**FD**) schemes to solve the *diffusive part* of the time-dependent or steady-state *temperature equation* including radioactive heating, in [1-D](./man/DiffOneD.md) and [2-D](./man/DiffTwoD.md). The solvers are located in [src/HeatEquation](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/HeatEquation/). So far, only *Dirichlet* and *Neumann* thermal boundary conditions are available. Most of the functions assume constant thermal parameters (except for the 1-D solvers and the 2-D defect correction solver). 

### [Heat Advection Equation](./man/AdvectMain.md)

The ```GeoModBox.jl``` provides different methods to advect certain properties within the model domain. The corresponding routines are structured in such a way, that any property can be advected with the described advection solvers, as long as the property is defined on the *centroids* including *ghost nodes* at all boundaries. Using passive tracers, one can, so far, choose to either advect the absolute temperature or the phase ID. 

## [Momentum Conservation Equation](./man/MomentumMain.md)

## Code Structure

### Initial Conditions 

### Scaling

## [Benchmarks and Examples](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/)

### [Gaussian Temperature Diffusion](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/DiffusionEquation/2D/Gaussian_Diffusion.jl)

![GaussianDiffusion](./assets/Gaussian_Diffusion_CNA_nx_100_ny_100.gif)

**Figure 1. Gaussian Diffusion.** Time-dependent, diffusive solution of a 2-D Gaussian temperature anomaly using the [Crank-Nicholson approach](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/HeatEquation/CNA.jl) in comparison to its analytical solution. Top Left: 2-D temperature field of the numerical solution and isotherms lines of the numerical (solid black) and analytical (dashed yellow) solution. Top Right: Total deviation to the analytical solution. Bottom Left: 1-D y-profile along x=0. Bottom Right: Root Mean Square total devation of the temperature over time. 

![GDResTest](./assets/Gaussian_ResTest.png)

**Figure 2. Resolution test.** Maximum *RMS* $\varepsilon$, maximum, and mean temperature for each **FD**-scheme and multiple resolutions. 

### [Rigid-Body-Rotation](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/AdvectionEquation/2D_Advection.jl)

![RigidBodyI](./assets/2D_advection_circle_RigidBody_upwind_100_100_nth_1.gif)

![RigidBodyII](./assets/2D_advection_circle_RigidBody_semilag_100_100_nth_1.gif)

![RigidBodyIII](./assets/2D_advection_circle_RigidBody_tracers_100_100_nth_1.gif)

**Figure 3. Rigid-Body-Rotation.** Time-dependent solution of a rotating circular temperature anomaly using the **upwind (first)**, **semi-lagrangian (second)**, and **tracer (third)** method. Within a circular area of our model domain the velocity is set to the velocity of a rigid rotation and outside euqal to zero. The temperature is scaled by the maximum temperature of the anomaly. 

### [Falling Block](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/StokesEquation/2D/FallingBlockBenchmark.jl)

![FallingBlockTD](./assets/Falling_block_ηr_0.0_tracers.gif)

**Figure 4. Isoviscous Falling Block.** Time-dependent solution of an isoviscous falling block example. The problem is solved with a solver for variable viscosities. The tracers advect the phase ID, which is used to interpolate the density and viscosity on the centroids and vertices, respectively. 

![FBSinkinVeloc](./assets/FallingBlock_SinkingVeloc_tracers.png)

**Figure 5. Falling Block Sinking Velocity.** Sinking velocity of the block with respect to the viscosity ratio $\eta_r$ at the initial condition. 

![FBFinalStage](./assets/FallingBlock_FinalStage_tracers.png)

**Figure 6. Falling Block Benchmark.** Final tracers distribution for specific cases with $\eta_r \ge 0 $. 

------------------
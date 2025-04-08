[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://geosci-ffm.github.io/GeoModBox.jl/)

# GeoModBox.jl
The **Geod**ynamic **Mod**elling Tool**Box** is a julia package mainly used for teaching purposes. The package provides different finite difference, staggered, discretization schemes to numerically solve the governing equations for a two-dimensional geodynamic problem. The governing equations are the conservation equations of 

1) [**energy**](./docs/src/man/DiffMain.md), 
2) [**momentum**](./docs/src/man/MomentumMain.md), 
3)  [**mass** and **compositon**](./docs/src/man/AdvectMain.md). 

The ```GeoModBox.jl``` includes a series of [exercises](/exercises/) and [examples](/examples/) of different geodynamically well defined problems. The exercises are given as Jupyter notebooks for the students to complete. The theoretical background and detailed explenations of the examples and functions are mainly given in the [documentation](https://lukasfuchs.github.io/GeoModBox.jl/).

<!-- 
- Different properties can be advected! 
### Discretization Schemes
#### Upwind 
#### Staggered Leap Frog (SLF)
#### Semi-Lagrangian
#### Passive Tracers -->

<!-- ## [Momentum Conservation Equation](./examples/StokesEquation/README.md) -->

<!--
- Constant viscosity
- Variable viscosity
- Coefficients assembly 
- Right-hand side updates
- Direct solution 
- Defect correction
-->

## Structure

<!-- The ```GeoModBox.jl``` is  -->

<!-- ## Initial Conditions  -->

<!---
- Set Up, Geometry, Grid, etc. 
- IniTemperature, IniVelocity, IniPhase,
- Tracer Initialization
- Solution of Stokes equation, Diffusion equation, Advection Equation, 
- ...
-->

<!-- ## Scaling -->


## [Benchmarks and Examples](./examples/)

### [Gaussian Temperature Diffusion](./examples/DiffusionEquation/2D/Gaussian_Diffusion.jl)
<img src="./examples/DiffusionEquation/2D/Results/Gaussian_Diffusion_CNA_nx_100_ny_100.gif" alt="drawing" width="600"/> <br>
**Figure 1. Gaussian Diffusion.** Time-dependent, diffusive solution of a 2-D Gaussian temperature anomaly using the [Crank-Nicholson approach](./src/HeatEquation/CNA.jl) in comparison to its analytical solution. Top Left: 2-D temperature field of the numerical solution and isotherms lines of the numerical (solid black) and analytical (dashed yellow) solution. Top Right: Total deviation to the analytical solution. Bottom Left: 1-D y-profile along x=0. Bottom Right: Root Mean Square total devation of the temperature over time. 

<img src="./examples/DiffusionEquation/2D/Results/Gaussian_ResTest.png" alt="drawing" width="600"/> <br>
**Figure 2. Resolution test.** Maximum *RMS* $\varepsilon$, maximum, and mean temperature for each **FD**-scheme and multiple resolutions. 

### [Rigid-Body-Rotation](./examples/AdvectionEquation/2D_Advection.jl)

<img src="./examples/AdvectionEquation/Results/2D_advection_circle_RigidBody_upwind_100_100_nth_1.gif" alt="drawing" width="300"/> 
<img src="./examples/AdvectionEquation/Results/2D_advection_circle_RigidBody_semilag_100_100_nth_1.gif" alt="drawing" width="300"/> 
<img src="./examples/AdvectionEquation/Results/2D_advection_circle_RigidBody_tracers_100_100_nth_1.gif" alt="drawing" width="300"/> 
<br>

**Figure 3. Rigid-Body-Rotation.** Time-dependent solution of a rotating circular temperature anomaly using the **upwind (first)**, **semi-lagrangian (second)**, and **tracer (third)** method. Within a circular area of our model domain the velocity is set to the velocity of a rigid rotation and outside euqal to zero. The temperature is scaled by the maximum temperature of the anomaly. 

### [Falling Block](./examples/StokesEquation/2D/FallingBlockBenchmark.jl)

<img src="./examples/StokesEquation/2D/Results/Falling_block_ηr_0.0_tracers.gif" alt="drawing" width="500"/> <br>

**Figure 4. Isoviscous Falling Block.** Time-dependent solution of an isoviscous falling block example. The problem is solved with a solver for variable viscosities. The tracers advect the phase ID, which is used to interpolate the density and viscosity on the centroids and vertices, respectively. 

<img src="./examples/StokesEquation/2D/Results/FallingBlock_SinkingVeloc_tracers.png" alt="drawing" width="500"/> <br>

**Figure 5. Falling Block Sinking Velocity.** Sinking velocity of the block with respect to the viscosity ratio $\eta_r$ at the initial condition. 

<img src="./examples/StokesEquation/2D/Results/FallingBlock_FinalStage_tracers.png" alt="drawing" width="500"/> <br>

**Figure 6. Falling Block Benchmark.** Final tracers distribution for specific cases with $\eta_r \ge 0$. 

<!--
- Blanckenbach
- Channel Flow
- Falling Block
- Gauss Diffusion, check! 
- RTI 
- Rigid Body Rotation, check! 
- Viscous Inclusion
-->

# Examples and Benchmarks

The `GeoModBox.jl` provides various one- and two-dimensional examples and benchmark problems for each of the governing equations. The examples demonstrate how to implement different numerical solvers, apply scaling, and evaluate the advantages and limitations of various finite difference schemes.

By clicking on the title of each document page, you will be directed to the corresponding Julia script in the [examples directory](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples).

**Advection** 
- [2-D advection with constant velocity field](./Advection/Advection2D.md)
- [Resolution test of 2-D advection](./Advection/AdvectionRestest2D.md)

**Heat Diffusion** 
- [1-D continental geotherm](./Diffusion/ContinentalGeotherm.md)
- [Comparison of FD schemes on a Gaussian anomaly](./Diffusion/GaussianDiffusion1D.md)
- [1-D oceanic geotherm](./Diffusion/OceanicGeotherm.md)
- [2-D Gaussian anomaly using iterative backward Euler solver](./Diffusion/BackwardEuler_DC.md)
- [2-D Gaussian anomaly using iterative forward Euler solver](./Diffusion/ForwardEuler_DC.md)
- [2-D resolution test with Gaussian anomaly](./Diffusion/GaussianDiffusion2D.md)
- [2-D Poisson equation resolution test](./Diffusion/PoissonRestest.md)
- [2-D Poisson equation with variable thermal properties](./Diffusion/PoissonVariablek.md)

**Thermal Convection Models** 
- [Bottom heated, isoviscous thermal convection](./BottomHeatedConvection.md)
- [Internally heated, isoviscous thermal convection](./Convection/InternallyHeatedConvection.md)
- [Mixed heated, isoviscous thermal convection](./Convection/MixedHeatedConvection.md)

**Stokes Equation**
- [1D channel flow with constant and depth-dependent viscosity](./Stokes/ChannelFlow1D.md)
- [2D falling block benchmark (direct or defect correction)](./Stokes/FallingBlockBenchmark.md)
- [2D falling block with constant or variable viscosity (defect correction)](./Stokes/FallingBlockDC.md)
- [2D Rayleigh–Taylor instability](./Stokes/RTI.md)
- [2D Rayleigh–Taylor instability benchmark and resolution test](./Stokes/RTI_growth_rate.md)
- [2D viscous inclusion problem](./Stokes/ViscousInclusion.md)
- [Van Keken Benchmark](./Stokes/VanKekenBenchmark.md)

**Thermo-Mechanical Shear Localization**
- Viscous shear localization

In the following, the runtime for each of the provided examples is listed as a reference. 

| Example                           | Total Runtime                                         |
| :-------------------------------- | :---------------------------------------------------- |
| **=== Advection ===**                                                                     |
| 2D_Advection.jl                   | 1) Upwind: 173 s                                      |
|    Resolution: 100x100            | 2) SLF: 173 s                                         |
|                                   | 3) Semi-lag: 177 s                                    |
|                                   | 4) Tracers: 459 s                                     |
| 2D_Advection_ResolutionTest.jl    | 1) save_fig = 1: 1.04 h                               |
|                                   | 2) save_fig = -1: 1109 s                              |
| **=== Heat Diffusion ===**                                                                |
| *1D* ---                                                                                  |
| ContinentalGeotherm_1D.jl         | 7.22 s                                                | 
| Heat_1D_discretization.jl         | 3.66 s                                                |
| OceanicGeotherm_1D.jl             | 691 ms                                                |
| *2D* ---                                                                                  |
| BackwardEuler.jl                  | 14.4 s                                                | 
| ForwardEuler.jl                   | 5.08 s                                                |
| Gaussian_Diffusion.jl             | 422 s                                                 | 
| Poisson_RestTest.jl               | 23.5 s                                                |
| Poisson_variable_k.jl             | 7.69 s                                                |
| **=== Stokes Equation ===**                                                               |
| *1D* ---                                                                                  |
| ChannelFlow_1D.jl                 | 1.53 s                                                |
| *2D* ---                                                                                  |
| FallingBlockBenchmark_*.jl        | 1) Steady State: 6.33 s (dc: 5.98 s)                  |
|                                   | 2) Time-Dependent                                     |
|                                   |     a) Upwind: 114 s (dc: 114 s)                      |
|                                   |     b) SLF: 68.3 s (dc: 68.8 s)                       |
|                                   |     c) Semi-lag: 319 s (dc: 325 s)                    |
|                                   |     d) Tracers: 326 s (dc: 328 s)                     |
| FallingBlockConstEta_DC.jl        | 379 ms                                                | 
| FallingBlockVarEta_DC.jl          | 40.3 s                                                |
| RTI.jl                            | 117 s                                                 |
| RTI_GrowthRate.jl                 | 46.6 s                                                |
| RTI_Growth_Rate_Res_Test_CNC.jl   | 414 s                                                 |
| RTI_Growth_Rate_Res_Test_CNM.jl   | 1573 s                                                |
| ViscousInclusion.jl               | 810 ms                                                |
| **=== Thermal Convection Models ===**                                                     |
| BottomHeated.jl                   |                                                 |
| InternallyHeated.jl               |                                                 |
| MixedHeated.jl                    |                                                       |

> **Note:** In `GeoModBox.jl`, thermal and kinematic boundary conditions are explicitly implemented within the solvers. The absolute values at *ghost nodes* are computed based on the values provided in the boundary condition tuple `BC`. Each tuple specifies the `type` (either Dirichlet or Neumann) and the corresponding `val`ue at each boundary.  
> For constant velocity boundary conditions, additional values must be defined in `val` for the boundary nodes (e.g., `BC.val.vxW`, `BC.val.vxE`). These additional values are required to **directly** solve the momentum equation using a the backslash operator and to update the right-hand side of the linear system. Furthermore, if non-zero, these values must be assigned to the initial boundary nodes of the respective velocity fields.  
> For more details on the implementation of constant velocity boundaries, refer to the documentation of the [viscous inclusion](./Stokes/ViscousInclusion.md) example or the [initial velocity setup](../Ini.md).

> **Note:** By default, the results of time-dependent examples in `GeoModBox.jl` are stored as GIF animations. To visualize solutions at specific time steps without generating a GIF, set the parameter `save_fig = 0`. In this case, individual plots are not saved, so caution is advised when running problems that require multiple time step iterations.

> **Note:** Some examples use *named tuples* to define constants and parameters. Alternatively, *mutable structures* can be used—particularly useful when parameters need to be modified after initialization (e.g., for scaling purposes). A full transition from *named tuples* to *mutable structures* is planned for future versions of `GeoModBox.jl`.
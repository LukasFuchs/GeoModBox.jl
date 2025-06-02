# Examples and Benchmarks

In the following, one can find different examples and benchmarks in 1-D and 2-D for each of the governing equations. The examples highlight how to implement the different solvers, how to utilize the scaling, and advantages and disadvantages of each finite difference scheme. 

By clicking on the title of each page one is directly directed to the julia example file within the [example directory](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples)

> **Note:** Within ```GeoModBox.jl``` the thermal and kinematic boundary conditions are explicitly implemented in the solvers, where the absolute values for the *ghost nodes* are calculated depending on the values given in the tuple ```BC```. Within the tuple, we define the ```type``` (Dirichlet or Neumann) and the corresponding ```val```ue at each boundary. 

> **Note:** Usually, the results of each example within the ```GeoModBox.jl``` are stored in *gif* animations (if it is a time-dependent problem). If one wants to plot the solution for certain time steps the parameter ```save_fig``` needs to be set to 0. This setting does not result in the generation of a gif file and the single plots are not saved! Thus, care needs to be taken if the problem needs multiple time step iterations. 

> **Note:** Some examples use *named tuples* to define the different constants and variables. Alternatively, *mutable structures* can also be used to define those parameters in the ```GeoModBox.jl```. The mutable structure are especially benificial, if one needs to edit the parameters after the have been defined, e.g., for scaling. 

## 2D

### [Backward Euler](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/HeatEquation/2D/BackwardEuler.jl)
-> A gaussian diffusion using the defection correction method. The results are compared to the analytical solution.  

### [Forward Euler](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/HeatEquation/2D/ForwardEuler.jl)
-> A gaussian diffusion using an explicit formulation. The results are compared to the analytical solution. 

### [Gaussian Diffusion](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/HeatEquation/2D/Gaussian_Diffusion.jl)
-> Script to solve a 2-D gaussian diffusion using multiple different finite difference schemes and resolutions. The results are compared to the analytical solution and the final plot shows a resolution test including each finite difference scheme.

### [Resolution Test Poisson Problem](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/HeatEquation/2D/Poisson_ResTest.jl)
-> Resolution test for a 2-D Poisson Problem. 

### [Poisson Variable Parameters](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/HeatEquation/2D/Poisson_variable_k.jl) 
-> 2-D Poisson Problem with variable thermal parameters.

# Advection Equation

## [2D Advection of an anomaly](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/AdvectionEquation/2D_Advection.jl)
-> ...

### [Resolution test of the 2-D advection](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/AdvectionEquation/2D_Advection_ResolutionTest.jl)
-> ...

# Momentum Equation

## 1D

### [Channel Flow (DC)](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/StokesEquation/1D/ChannelFlow_1D.jl)
<!--
### Channel Flow
&emsp; Assuming the horizontal pressure gradient is constant and flow within a channel is only driven by the pressure and/or by a constant horizonal velocity at the surface (or at the bottom, or both), the stokes equation describes the horizontal flow velocity within the channel and simplifies to: 

For the given setup I can assume that the vertical velocity is zero and thus equation (3) simplifies to the last expression of equation (2).

&emsp; This directory contains a script to calculate the horizontal velocity for a two-dimensional Couette(-Poiseuille) channel flow with constant and logarithmically, with depth varying viscosity and to compare the numerical solution with its analytical solution. The depth-dependent viscosity is defined as: 

$\eta = \eta_0 exp(log(m) \frac{H-z}{H})$,

where *m* is the viscosity ratio of $\frac{\eta_1}{\eta_0}$, $\eta_0$ and $\eta_1$ are the bottom and surface viscosities, respectively, *H* is the model height, and *z* the depth. 

&emsp;Considering the definition of the viscosity as given in equation (4), one can derive an analytical solution of the horizontal velocity from the 1-D stokes equation in *x*-direction by twice integrating equation (1). The analytical solution with depth depends on the viscosity ratio *m*, the horizontal pressure gradient $\frac{\partial P}{\partial x}$, and the shear velocity at the surface $v_{x,0}$. For an upward pointing coordinate system (*z* positive) the analytical solution is given as: 

$v_{x,ana} =-\frac{1}{2 \eta_0} \frac{\partial P}{x} (Hz - z^2) + v_{x,0}\frac{z}{H}$,&emsp;&emsp; if $m = 1$, and &emsp;&emsp;&emsp; (5)

$v_{x,ana} = -\frac{\partial P}{\partial x} \frac{H}{\eta_0 log(m)} (\frac{m^{-\frac{z}{H}}}{m-1}(z(m-1)+H) - \frac{H}{m-1})-m^{-\frac{z}{H}} m \frac{v_{x,0}}{m-1} + \frac{v_{x,0}m}{m-1}$, &emsp;&emsp; if $m \neq 0$.&emsp;&emsp;&emsp; (6)

&emsp;The numerical solution is calculated using fixed boundary velocities, which are defined by the analytical solution of the horizontal velocity as defined in equations (5) and (6) and I simply flip the analytical solution so that it fits to the downward point coordinate system I use in the code. -->

## 2D

### [FallingBlockBenchmark()](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/StokesEquation/2D/FallingBlockBenchmark.jl)
-> A sript, solving the falling block benchmark for a viscosity range from -6 to 6 order of magnitude. The script stores the [sinking velocity](../assets/FallingBlock_SinkingVeloc_tracers.png) of the block at the initial configuration and the [final marker distribution](../assets/FallingBlock_FinalStage_tracers.png) for models with a viscosity ration over and equal 0. 

### [FallingBlockConstEta_DC()](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/StokesEquation/2D/FallingBlockConstEta_Dc.jl) 
-> A sript, solving the falling block problem assuming a constant viscosity and using the defect correction method.  

### [FallingBlockVarEta_DC()](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/StokesEquation/2D/FallingBlockVarEta_DC.jl)
-> A sript, solving the falling block problem assuming a variable viscosity and using the defect correction method. The advection is only provided with tracers

### [RTI](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/StokesEquation/2D/RTI.jl)
-> To be added (*tba*)

### [Viscous Inclusion](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/StokesEquation/2D/ViscousInclusion.jl)
-> To be added (*tba*)

# Thermal 
-> To be added (*tba*)

## [Bottom Heated Convection]()
-> To be added (*tba*)

## [Internally Heated Convection]()
-> To be added (*tba*)

## [Mixed Heated Convection]()
-> To be added (*tba*)




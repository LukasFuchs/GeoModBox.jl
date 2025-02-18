# GeoModBox.jl
**Geod**ynamic **Mod**elling Tool**Box** is a julia package mainly used for teaching purposes. The package provides different finite differences, fully staggered, discretization schemes to numerically solve the governing equations for a two-dimensional geodynamic problem. The governing equations are the conservation equations of 
1) [**energy**](/examples/HeatEquation/README.md), 
2) [**momentum**](/examples/...), 
3)  [**mass** and **compositon**](/examples/AdvectionEquation/README.md). 

&emsp; The ```GeoModBox.jl``` includes a series of [exercises](./exercises/) and [examples](./examples/) of different geodynamically well defined problems. The exercises are given as Jupyter notebooks for the students to complete. Most of the theoretical background is given in the corresponding README.md files and in the documentation. 

<!-- I will add some question directly here for the time beeing: I would like to add similar informations to this package as I did in my [FDCSGm repository](https://github.com/LukasFuchs/FDCSGm). However, I am not familiar with the doc option in github yet, that's why I added all the information in a README.md file for each module etc. I guess the doc option would be more suitable to add the details and simply add some figures and general information in the README.md. Have you ever worked with doc in github? -->

------------------
------------------
## Staggered Finite Difference 
<!--- brief general description of a staggered finite difference scheme
- Energy equation
- Momentum equation-->
------------------
------------------
## Energy equation 
<!-- - General formulation of the energy equation (advection + diffusion term + radiogenic heating)
- Operator splitting
- Solution of the conductive part -->
------------------
### Heat Diffusion
<!--- 
- Conducitve part of the energy equation + radiogenic heating + ?
- Constant thermal parameters
- Variable thermal paramters 
------------------
<!-- #### Numerical Schemes


All numerical schemes methods can be used in the [thermal convection code]() and the [Blankenbach Benchmark]() and are generally available to chose in the code. A more detailed analysis on the accuracy of each discretization scheme and the effect of the grid resolution is given in the [Gaussian Diffusion Benchmark](). 

##### Explicit, FTCS
##### Implicit, FTCS
##### Cranck-Nicolson Approach (CNA)
##### Alternating Direction Implicit (ADI) 
------------------
##### Thermal Boundary Conditions

1. Dirichlet
2. Neumann
------------------
#### Steady State Solution
##### Poisson solution, constant *k*
##### Poisson solution, variable *k* -->
### Heat Advection 
<!-- 
### Discretization Schemes
#### Upwind 
#### Staggered Leap Frog (SLF)
#### Semi-Lagrangian
#### Passive Tracers -->
------------------
------------------
## Momentum Equation
------------------
------------------
## Benchmarks
------------------
------------------
## Examples
------------------
------------------

# GeoModBox.jl
**Geodynamic Modelling** Tool**Box** is a julia package mainly used for teaching purposes. The package provides different finite differences, fully staggered, discretization schemes for the governing equations to describes a geodynamic problem. The geoverning equations are the conservation equations of **energy**, **mass**, **momentum**, and **compositon**. 

I will add some question directly here for the time beeing: I would like to add similar informations to this package as I did in my [FDCSGm repository](https://github.com/LukasFuchs/FDCSGm). However, I am not familiar with the doc option in github yet, that's why I added all the information in a README.md file for each module etc. I guess the doc option would be more suitable to add the details and simply add some figures and general information in the README.md. Have you ever worked with doc in github? 

------------------
## Energy equation 
### Heat Diffusion
------------------
#### Discretization Methods
##### Explicit, FTCS
##### Implicit, FTCS
##### Cranck-Nicolson Approach (CNA)
##### Alternating Direction Implicit (ADI)
------------------
#### Steady State Solution
##### Poisson solution, constant *k*
##### Poisson solution, variable *k*
------------------
------------------
## Advection Equation
### Discretization Schemes
#### Upwind 
#### Staggered Leap Frog (SLF)
#### Semi-Lagrangian
#### Passive Tracers
------------------
------------------
## Stokes Equation
------------------
------------------
## Bechmarks
------------------
------------------
## Examples
------------------
------------------

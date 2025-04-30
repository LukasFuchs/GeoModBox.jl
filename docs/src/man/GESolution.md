# Solving Differential Equations

- different methods to solve PDE -> Finit differences
- different methods to discretize equations 
- different methods to solve discretized equations

## Staggered Finite Difference

- conservation of stresses, heat-flux, 
- required for variable thermal parameters
- Total grid 

# Thermal convection

The given equations can be used to solve for the pressure and velocity of a 2-D thermal convection (awaiting further minor implementations also for variable thermodynamic parameters). However, it is quite common to simplify the governing euqations using certain approximations. 

## Approximations 

A very prominent approximation for a thermal convection is the so called *Boussinesq* approximation. Thereby, one assumes that all thermodynamic parameters are constant and one neglects adiabatic temperature effects. ... buoyancy term ...

### Equation of State

In case of a thermal convection the density is considered to be temperature dependent. Thus, the buoyance term on the right-hand side of equation $(7)$ is temperature dependent (and pressure, but one can neglect this effect here so far) and can be approximated with the so-called *equation of state* for the density. Here, its is a linear approximation of the change of density due to temperature variations and can be defined as:

$\begin{equation}
\rho = \rho_0 (1-\alpha T),
\end{equation}$

where $ρ_0$ is the reference density and $\alpha$ the thermal expansion coefficient [ $1/K$ ]. 

- change of equations (boussinesq approximation)

# Scaling

- scaling constants
- non-dimensional equations

... *tba* ...

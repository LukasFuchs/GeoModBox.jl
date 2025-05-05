# Solving Differential Equations

- different methods to solve PDE -> Finit differences
- different methods to discretize equations 
- different methods to solve discretized equations

## Staggered Finite Difference

- conservation of stresses, heat-flux, 
- required for variable thermal parameters
- Total grid 

# Thermal convection

The given equations can be used to solve for the pressure and velocity of a 2-D thermal convection (awaiting further minor implementations also for variable thermodynamic parameters, $\rho$, $c_p$, and $k$). However, it is quite common to simplify the governing euqations using certain approximations. 

## Approximations 

A very prominent approximation for a thermal convection is the so called *Boussinesq* approximation. Thereby, one assumes that all thermodynamic parameters are constant and one neglects adiabatic temperature effects within the *temperature equation*. The spatial density variations are assumed to be very small and are only taken into account in the buoyancy term of the *momentum equation*. Thereby, one assumes the density is only dependent on the temperature and can be defined by an *equation of state*.

### Equation of State

One can define different versions of the *equation of state*, where in this case, its is assumed to be a linear approximation of the change of density due to temperature and can be defined as:

$\begin{equation}
\rho = \rho_0 (1-\alpha T),
\end{equation}$

where $ρ_0$ is the reference density and $\alpha$ the thermal expansion coefficient [ $1/K$ ]. 

... definition total pressure ...

The corresponding equations to solve a thermal convection using a *Boussinesq* approsimation are then: 

**Momentum equation**

*x-component*

$\begin{equation}
0 = 
\end{equation}$

*y-component*

$\begin{equation}

\end{equation}$

**Temperature equation**

$\begin{equation}

\end{equation}$

- change of equations (boussinesq approximation)

The *continuity equation* remains the same. 

*...example...*

# Scaling

$\begin{equation}\begin{split}

\end{split}\end{equation}$

- scaling constants
- non-dimensional equations

... *tba* ...

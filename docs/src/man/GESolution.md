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

One can define different versions of the *equation of state*. For the here analysed thermal convection, its is a linear approximation of the change of density due to temperature and defined as:

$\begin{equation}
\rho = \rho_0 (1-\alpha T),
\end{equation}$

where $ρ_0$ is the reference density and $\alpha$ the thermal expansion coefficient [ $1/K$ ]. 

Replacing the buoyancy term in the right hand side of the *momentum equation* and considering the definition of the total pressure of

$\begin{equation}
P_t = P_{dyn} + P_{hydr},
\end{equation}$

where $P_{dyn}$ and $P_{hydr}$ are the dynamic and hydrostatic pressure, respectively,and the definition of the hydrostatic pressure gradient of

$\begin{equation}
\frac{\partial{P_{hydr}}}{\partial{y}}=\rho_0 g,
\end{equation}$

results in a different form of the $y$-component of the dimensional *momentum equation*.

**Governing equations**

The corresponding dimensional equations to solve a thermal convection using a *Boussinesq* approsimation are: 

**Momentum equation**

*x-component*

$\begin{equation}
0 = -\frac{\partial{P_{dyn}}}{\partial{x}}+\frac{\partial{\tau_{xj}}}{\partial{x_j}},
\end{equation}$

*y-component*

$\begin{equation}
0 = -\frac{\partial{P_{dyn}}}{\partial{y}}+\frac{\partial{\tau_{yj}}}{\partial{x_j}} - \rho_0 g \alpha T,
\end{equation}$

where $P_{dyn}$ is the dynamic pressure in [ $Pa$ ], $\tau_{ij}$ is the deviatoric stress tensor in [ Pa ], $\rho_0$ is the reference density in [ $kg/m^3$ ], $g$ is the gravitational acceleration in [ $m/s^2$ ], $\alpha$ is the thermal expansion coefficient in [ $1/K$ ], and $T$ the absolute temperature in [ $K$ ]. 

**Temperature equation**

$\begin{equation}
\left(\frac{\partial{T}}{\partial{t}} + v_j \frac{\partial{T}}{\partial{x_j}}\right) = - \kappa \frac{\partial^2{T}}{\partial{x^2_i}} + \frac{Q}{\rho_0 c_p},
\end{equation}$

where $t$ is the time in [ $s$ ], $v_j$ is the velocity in the $j$-th direction in [ $m/s$ ], $\kappa = k/\rho/c_p$ is the thermal diffusivity in [ $m^2/s$ ], $Q$ is the heat production rate per volume in [ $W/m^3$ ], and $c_p$ is the specific heat capacity in [ $J/kg/K$ ]. To see how this is implemented wiwthin the code, please see the [thermal convection example codes](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/MixedHeatedConvection/).

**Continuity equation**

$\begin{equation}
\frac{\partial{v_i}}{\partial{x_i}} = 0.
\end{equation}$

# Scaling



$\begin{equation}\begin{split}

\end{split}\end{equation}$

- scaling constants
- non-dimensional equations

... *tba* ...

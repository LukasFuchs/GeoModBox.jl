# Solving Differential Equations

## Governing Equations

The general governing equations for solving a geodynamical problem, neglecting adiabatic effects and assuming only radioactive heat source, are the conservation equations of 

**Momentum**

$\begin{equation}
\rho \left(\frac{\partial{v_{i}}}{\partial{t}} + v_{j}\frac{\partial{v_{i}}}{\partial{x_{j}}}\right) = -\frac{\partial{P}}{\partial{x_{i}}} + \frac{\tau_{ij}}{\partial{x_j}} + \rho g_{i},
\end{equation}$

where 
$\rho$ is the density [kg/m³], 
$v_i$ is the velocity in the $i$-th direction [m/s],
$t$ is time [s],
$\partial/\partial{t}$ is the time derivative, 
$\partial/\partial x_i$ is a directional derivative in $i$, 
$P$ the total pressure [Pa], 
$\tau_{ij}$ is the deviatoric stress tensor [Pa], and 
$\boldsymbol{g}$ is the gravitational acceleration vector [m/s²]. 

**Energy**

$\begin{equation}
\rho c_p \left(\frac{\partial T}{\partial t} + v_j\frac{\partial{T}}{\partial{x_j}}\right) = -\frac{\partial q_i}{\partial x_i} + \rho H,
\end{equation}$

where 
$c_p$ is the specific heat capacity [J/kg/K],
$T$ is temperature [K],
$q_i$ is the heat flux [W/m²] in direction $i$,
$H$ is the internal heat production per unit mass [W/kg]. 

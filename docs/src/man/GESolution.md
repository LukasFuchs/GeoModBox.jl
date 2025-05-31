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

**Mass**

$\begin{equation}
\frac{\partial{v_i}}{\partial{x_i}} = 0.
\end{equation}$

Repeated indices imply summation.

---

Ordinary and partial differential equations (ODEs and PDEs) can be solved through various approaches—occasionally *analytically*, but more commonly *numerically* due to their inherent complexity. Among numerical methods, prominent techniques include *integral*-based methods, such as the *finite element* and *spectral* methods, as well as the *finite difference* method.

The ```GeoModBox.jl``` framework primarily employs the **finite difference method**. While each numerical approach has its own strengths and limitations, the choice often depends on the user's familiarity and comfort with the method. Nonetheless, the finite difference method is relatively straightforward and pedagogically advantageous, as its discretized form closely resembles the original differential equations. Furthermore, it is computationally efficient, making it well-suited for performance-critical applications.

In general, the finite difference method aims to approximate **differential operators** using finite differences derived from a Taylor series expansion. For further details, refer to the [lecture notes](https://lukasfuchs.wordpress.com/numerical-methods-in-geophysics/) or see the reference below.

## Staggered Finite Difference

To solve differential equations within a given domain using the finite difference method, it is first necessary to generate a *numerical grid* on which finite differences can be computed. The most straightforward approach is to discretize the domain using a *regular*, *uniform* grid, where the spacing between grid points is constant and all variables are defined at the same locations. Such grids are commonly used to solve equations like the Poisson equation, the heat equation, or advective transport equations.

However, in many cases, certain limitations or physical constraints require a different arrangement of variable locations to ensure numerical stability and the conservation of physical properties. For instance, solving the *momentum equation* with **variable viscosity** typically requires a *fully staggered grid* to correctly preserve stress continuity across adjacent grid points. A similar consideration applies to the *temperature equation* when using **variable thermal conductivity**.

Staggered grids also offer advantages in implementing boundary conditions. For example, with *Neumann* thermal boundary conditions, the heat flux across a boundary can be naturally evaluated at staggered points using so-called *ghost points*. These ghost points can also be employed to impose *Dirichlet* boundary conditions. This approach helps maintain consistent accuracy and order of the finite difference scheme both at internal and boundary points. 

For these reasons, `GeoModBox.jl` adopts a staggered grid for solving the *temperature equation*. The complete grid structure used for the governing equations in `GeoModBox.jl` is illustrated below:

![2DGrid_total](../assets/Grid_total.png)

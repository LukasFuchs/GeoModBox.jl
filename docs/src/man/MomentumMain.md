# Momentum Equation

On geological time scales, Earth's mantle and lithosphere behave as a highly viscous fluid that moves and deforms in response to applied forces. In general, three major types of forces influence fluid motion: inertial, surface, and body (volumetric) forces. A commonly used expression describing fluid motion is given by:

$\begin{equation}
\rho \frac{D v_i}{Dt} = \frac{\partial{\sigma_{ij}}}{\partial{x_j}} + \rho g_i,
\end{equation}$

where  

$\rho$ is the density [kg/m³],  
$v_i$ is the velocity component [m/s] in direction $i$,  
$\sigma_{ij}$ is the Cauchy stress tensor [Pa],  
$g_i$ is the gravitational acceleration [m/s²] in direction $i$, and  
$\frac{D}{Dt}$ is the Lagrangian (material) time derivative, which in Eulerian form is written as:

$\begin{equation}
\frac{D}{Dt} = \frac{\partial{}}{\partial{t}} + v_j\frac{\partial}{\partial{x_j}}.
\end{equation}$

The Cauchy stress tensor is commonly decomposed into an isotropic pressure contribution and a deviatoric stress component:

$\begin{equation}
\sigma_{ij} = \tau_{ij} -P_t \delta_{ij},
\end{equation}$

where $P_t$ is the total pressure (including dynamic and hydrostatic components), $\tau_{ij}$ is the deviatoric stress tensor, and $\delta_{ij}$ is the Kronecker delta.

Substituting this decomposition into the momentum equation and expanding the material derivative yields the Navier–Stokes equation in Eulerian form:

$\begin{equation}
\rho \left(\frac{\partial{v_{i}}}{\partial{t}} + v_{j}\frac{\partial{v_{i}}}{\partial{x_{j}}}\right) = -\frac{\partial{P}}{\partial{x_{i}}} + \frac{\partial{\tau_{ij}}}{\partial{x_j}} + \rho g_{i},
\end{equation}$

where Einstein summation is implied over repeated indices.

## Constitutive Relation

To solve the momentum equation, the rheology of the material must be specified. For a purely viscous (Newtonian) fluid, the deviatoric stress tensor is related to the strain-rate tensor via:

$\begin{equation}
\tau_{ij} = 2 \eta \cdot \dot{\varepsilon}_{ij},
\end{equation}$

where  
$\eta$ is the dynamic viscosity [Pa·s] and  
$\dot{\varepsilon}_{ij}$ is the strain-rate tensor [1/s], defined as:

$\begin{equation}
\dot{\varepsilon}_{ij} = \frac{1}{2} \left(\frac{\partial{v_i}}{\partial{x_j}} + \frac{\partial{v_j}}{\partial{x_i}}\right).
\end{equation}$

## Stokes Equation

In geodynamic applications, inertial forces are typically negligible compared to viscous and gravitational forces. Under this Stokes flow approximation, the momentum equation simplifies to:

$\begin{equation}
0 = -\frac{\partial{P}}{\partial{x_{i}}} + \frac{\partial{\tau_{ij}}}{\partial{x_j}} + \rho g_{i},
\end{equation}$

which is commonly referred to as the Stokes equation. This approximation corresponds to flow at very low Reynolds numbers, typical for mantle convection.

Using the constitutive relation for a viscous fluid yields:

$\begin{equation}
0 = -\frac{\partial{P}}{\partial{x_{i}}} + \frac{\partial}{\partial{x_j}} \left( \eta \left(\frac{\partial{v_i}}{\partial{x_j}} + \frac{\partial{v_j}}{\partial{x_i}}\right)\right) + \rho g_{i}.
\end{equation}$

Assuming constant viscosity simplifies the equation further to:

$\begin{equation}
0 = -\frac{\partial{P}}{\partial{x_{i}}} + \eta \left(\frac{\partial^2{v_i}}{\partial{x_j^2}} + \frac{\partial^2{v_j}}{\partial{x_i^2}}\right) + \rho g_{i}.
\end{equation}$

In two dimensions, solving this system requires discretization of both the $x$- and $y$-components of the momentum equation, together with the mass conservation equation. For implementation details, refer to the [1D](./MomentumOneD.md) and [2D](./MomentumTwoD.md) momentum equation documentation.

## Continuity Equation

The Stokes system provides two equations for the three unknowns $v_x$, $v_y$, and $P$. To close the system, the continuity equation is used. Assuming an incompressible fluid (Boussinesq approximation), the mass conservation equation becomes:

$\begin{equation}
\frac{\partial{v_i}}{\partial{x_i}} = 0.
\end{equation}$

Together, the momentum and continuity equations form the Stokes system, which allows solving for the velocity components and pressure.

## Examples

The following examples demonstrate applications of the Stokes equations implemented in `GeoModBox.jl`:

- [1D channel flow with constant and depth-dependent viscosity](./examples/ChannelFlow1D.md)  
- [2D falling block benchmark](./examples/FallingBlockBenchmark.md)  
- [2D falling block with constant viscosity (defect correction)](./examples/FallingBlockDC.md)  
- [2D falling block with variable viscosity (defect correction)](./examples/FallingBlockDC.md)  
- [2D viscous inclusion problem](./examples/ViscousInclusion.md)  
- [2D Rayleigh–Taylor instability benchmark](./examples/RTI_growth_rate.md)

Examples of coupled temperature–momentum systems (i.e., thermal convection models) using operator splitting include:

- [Mixed heated convection models](./examples/MixedHeatedConvection.md)

## Exercises

- [Steady-state, isoviscous 2D falling block](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/09_2D_Falling_Block.ipynb)  
- [Time-dependent, isoviscous 2D falling block](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/10_2D_Falling_Block_td.ipynb)  
- [2D thermal convection (isoviscous)](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/11_2D_Thermal_Convection.ipynb)  
- [Scaled 2D thermal convection](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/12_2D_Thermal_Convection_scaled.ipynb)  
- [Blankenbach benchmark with resolution study](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/13_Blankenbach_Benchmark.ipynb)
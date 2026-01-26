# Momentum Equation

On geological time scales, Earth's mantle and lithosphere behave like a fluid that moves and deforms in response to forces. In general, three major types of forces influence fluid motion: **inertial**, **surface**, and **body (volumetric)** forces. A commonly used expression to describe such motion is:

$\begin{equation}
\rho \frac{D v_i}{Dt} = \frac{\partial{\sigma_{ij}}}{\partial{x_j}} + \rho g_i,
\end{equation}$

where 
$\rho$ is the density [kg/m³], 
$v_i$ is the velocity component [m/s] in direction $i$, 
$\sigma_{ij}$ is the **Cauchy stress tensor** [Pa], 
$g_i$ is the gravitational acceleration [m/s²] in direction $i$, and 
$\frac{D}{Dt}$ is the Lagrangian (material) time derivative, expressed in Eulerian form as: 

$\begin{equation}
\frac{D}{Dt} = \frac{\partial{}}{\partial{t}} + v_j\frac{\partial}{\partial{x_j}}.
\end{equation}$

The Cauchy stress tensor is commonly decomposed as:

$\begin{equation}
\sigma_{ij} = \tau_{ij} -P_t \delta_{ij},
\end{equation}$

where $P_t$ is the total pressure (including dynamic and hydrostatic components), $\tau_{ij}$ is the **deviatoric stress tensor**, and $\delta_{ij}$ is the Kronecker Delta.

Substituting into the momentum equation and expanding the material derivative yields the **Navier–Stokes equation** in Eulerian form:

$\begin{equation}
\rho \left(\frac{\partial{v_{i}}}{\partial{t}} + v_{j}\frac{\partial{v_{i}}}{\partial{x_{j}}}\right) = -\frac{\partial{P}}{\partial{x_{i}}} + \frac{\partial{\tau_{ij}}}{\partial{x_j}} + \rho g_{i},
\end{equation}$

with Einstein summation implied over repeated indices.

## Constitutive Relation

To solve the momentum equation, we must define the rheology of the material. For a **purely viscous** fluid, the deviatoric stress is related to the strain rate by:

$\begin{equation}
\tau_{ij} = 2 \eta \cdot \dot{\varepsilon}_{ij},
\end{equation}$

where 
$\eta$ is the dynamic viscosity [Pa·s] and 
$\dot{\varepsilon}_{ij}$ the **strain rate tensor** [1/s], given by: 

$\begin{equation}
\dot{\varepsilon}_{ij} = \frac{1}{2} \left(\frac{\partial{v_i}}{\partial{x_j}} + \frac{\partial{v_j}}{\partial{x_i}}\right).
\end{equation}$

## Stokes Equation

In the mantle, inertial forces are typically negligible compared to viscous and gravitational forces. Under this **Stokes flow** approximation, the momentum equation simplifies to:

$\begin{equation}
0 = -\frac{\partial{P}}{\partial{x_{i}}} + \frac{\partial{\tau_{ij}}}{\partial{x_j}} + \rho g_{i},
\end{equation}$

Using the constitutive relation, we obtain:

$\begin{equation}
0 = -\frac{\partial{P}}{\partial{x_{i}}} + \frac{\partial}{\partial{x_j}} \left( \eta \left(\frac{\partial{v_i}}{\partial{x_j}} + \frac{\partial{v_j}}{\partial{x_i}}\right)\right) + \rho g_{i}.
\end{equation}$

Assuming **constant viscosity** simplifies this further to:

$\begin{equation}
0 = -\frac{\partial{P}}{\partial{x_{i}}} + \eta \left(\frac{\partial^2{v_i}}{\partial{x_j^2}} + \frac{\partial^2{v_j}}{\partial{x_i^2}}\right) + \rho g_{i}.
\end{equation}$

Solving this system requires discretization of both the $x$- and $y$-components of the **momentum equation**, as well as the **mass conservation** equation. For implementation details, refer to the [1D](./MomentumOneD.md) and [2D](./MomentumTwoD.md) momentum equation documentation.

## Continuum Equation

The Stokes system provides two equations for the three unknowns $v_x$, $v_y$, and $P$. To close the system, we use the **continuity equation**. Assuming an incompressible fluid (Boussinesq approximation), the mass conservation equation becomes:

$\begin{equation}
\frac{\partial{v_i}}{\partial{x_i}} = 0.
\end{equation}$

Together, the momentum and continuity equations allow solving for $v_x$, $v_y$, and $P$.

## Examples

The following examples demonstrate the application of the Stokes equations:

- [1D channel flow with constant and depth-dependent viscosity](./examples/ChannelFlow1D.md)  
- [2D falling block benchmark](./examples/FallingBlockBenchmark.md)  
- [2D falling block with constant viscosity (defect correction)](./examples/FallingBlockDC.md)  
- [2D falling block with variable viscosity (defect correction)](./examples/FallingBlockDC.md)  
- [2D viscous inclusion problem](./examples/ViscousInclusion.md)  
- [2D Rayleigh–Taylor instability benchmark](./examples/RTI_growth_rate.md)

Examples of coupled temperature–momentum systems (i.e., **convection models**) using **operator splitting** include:

- [Mixed heated convection models](./examples/MixedHeatedConvection.md)


## Exercises

- [Steady-state, isoviscous 2D falling block](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/09_2D_Falling_Block.ipynb)  
- [Time-dependent, isoviscous 2D falling block](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/10_2D_Falling_Block_td.ipynb)
- [2D thermal convection (isoviscous)](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/11_2D_Thermal_Convection.ipynb)  
- [Scaled 2D thermal convection](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/12_2D_Thermal_Convection_scaled.ipynb)  
- [Blankenbach benchmark with resolution study](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/13_Blankenbach_Benchmark.ipynb)







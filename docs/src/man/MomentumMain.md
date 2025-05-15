# Momentum Equation

On geological time scales, Earth’s mantle and lithosphere do behave like a fluid and move and deform. A fluid does generally move due to forces acting on it whereas the forces must be in balance. In general, there are three major forces one might consider, i.e., inertia, surface, and volumetric forces. A common equation to describes such motion is given by: 

$\begin{equation}
\rho \frac{D \overrightharpoon{v}}{Dt} = \overrightharpoon{\nabla} \cdot \boldsymbol{\sigma} + \rho \boldsymbol{g},
\end{equation}$

where $\rho$ is the density [kg/m³], $\overrightharpoon{v}$ is the velocity vector [m/s], $\boldsymbol{\sigma}$ is the *Cauchy stress tensor* [Pa], $\boldsymbol{g}$ is the gravitational acceleration vector [m/s²], and the term on the left-hand side is the Lagrangian time derivative which is in Eulerian form $\frac{D}{Dt} = \frac{\partial{}}{\partial{t}} + \overrightharpoon{v} \cdot \overrightharpoon{\nabla}$. 

The *Cauchy stress tensor* is given by: 

$\begin{equation}
\boldsymbol{\sigma} = -\overrightharpoon{\nabla}{P} + \overrightharpoon{\nabla} \cdot \boldsymbol{\tau},
\end{equation}$

where $P$ is the total pressure ($P = P_{dynamic} + P_{hydrostatic}$) and $\boldsymbol{\tau}$ the *deviatoric stress tensor*. 

In Eulerian form, equation $(1)$ is given by (i.e., the **Navier-Stokes equation**):

$\begin{equation}
\rho \left(\frac{\partial{v_{i}}}{\partial{t}} + v_{j}\frac{v_{i}}{\partial{x_{j}}}\right) = -\frac{\partial{P}}{\partial{x_{i}}} + \frac{\tau_{ij}}{\partial{x_j}} + \rho g_{i},
\end{equation}$

where summation over repeated indices is implied.

# Constitutive Relation

To solve equation $(3)$, one needs to define a rheology which describes how a material deforms under certain applied forces. For a purely viscous medium, one can be define a constitutive relationship between stress and strain rate in the form of:

$\begin{equation}
\tau_{ij} = 2 \eta \cdot \dot{\varepsilon}_{ij},
\end{equation}$

where $\eta$ is the dynamic viscosity in [Pa s] and $\dot{\varepsilon}_{ij}$ the *strain rate tensor* in [1/s] and given by: 

$\begin{equation}
\dot{\varepsilon}_{ij} = \frac{1}{2} \left(\frac{\partial{v_i}}{\partial{x_j}} + \frac{\partial{v_j}}{\partial{x_i}}\right).
\end{equation}$

# Stokes Equation

Assuming that the inertia forces are negligible in comparison to the gravitational forces, which is generally the case for a high viscous medium like the Earth's mantle, one can further simplify equation $(3)$ to:

$\begin{equation}
0 = -\frac{\partial{P}}{\partial{x_{i}}} + \frac{\tau_{ij}}{\partial{x_j}} + \rho g_{i},
\end{equation}$

or in the form of the unknowns $v_x$, $v_y$, and $P$:

$\begin{equation}
0 = -\frac{\partial{P}}{\partial{x_{i}}} + \frac{\partial}{\partial{x_j}} \eta \left(\frac{\partial{v_i}}{\partial{x_j}} + \frac{\partial{v_j}}{\partial{x_i}}\right) + \rho g_{i}.
\end{equation}$

Assuming a constant viscosity further simplifies equation $(7)$ to: 

$\begin{equation}
0 = -\frac{\partial{P}}{\partial{x_{i}}} + \eta \left(\frac{\partial^2{v_i}}{\partial{x_j^2}} + \frac{\partial^2{v_j}}{\partial{x_i^2}}\right) + \rho g_{i}.
\end{equation}$

To numerically solve for the three unknows $v_x$, $v_y$, and $P$, one needs to discretize the $x$- and $y$*-component* of the **momentum equation** and the **continuum equation**. For more details on how this is done, please see the [1-D](./MomentumOneD.md) and [2-D](./MomentumTwoD.md) documentations of the momentum equation. 

# Continuum Equation

Equation $(8)$ provides us two equations for three unknowns. Thus, one also needs to consider the *mass conservation equation*, where one can further simplify the problem by assuming an incompressible fluid (i.e., Boussinesq-approximation):

$\begin{equation}
\frac{\partial{v_i}}{\partial{x_i}} = 0.
\end{equation}$

Equations $(8)$ and $(9)$ enable us to solve for the three unknowns $v_x$, $v_y$, and $P$. 

The examples for a pure *stokes equation* problem include: 
- [a 1-D channel flow problem for constant and depth-dependent viscosity](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/StokesEquation/1D/ChannelFlow_1D.jl)

- [a 2-D falling block benchmark](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/StokesEquation/2D/FallingBlockBenchmark.jl)

- [a 2-D falling block example with constant viscosity using the defect correction method](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/StokesEquation/2D/FallingBlockConstEta_DC.jl)

- [a 2-D falling block example with variable viscosity using the defect correction method](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/StokesEquation/2D/FallingBlockVarEta_DC.jl)

- [a 2-D viscous inclusion problem](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/StokesEquation/2D/ViscousInclusion.jl), and

- [a 2-D Rayleigh-Taylor instability benchmark](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/StokesEquation/2D/RTI.jl)

The exercises include: 

- [a steady state, isoviscous 2-D falling block](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/09_2D_Falling_Block.ipynb)

- [a time-dependent, isoviscous 2-D falling block](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/09_2D_Falling_Block_td.ipynb)

# Further applications

The examples include models of different [mixed heated convection system](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/MixedHeatedConvection/) to highlight the coupling between the *temperature* and *momentum equation* using the *operator splitting method*. Within these examples all methods to solve each equation can be applied. 

The exercises include an example to solve [an isovisous, 2-D thermal convection](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/11_2D_Thermal_Convection.ipynb) as well as [a scaled version thereof](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/12_2D_Thermal_Convection_scaled.ipynb). 

The final exercises fo complete the course is to complete the [Blankenbach benchmark including a resolution test](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/13_Blankenbach_Benchmark.ipynb). 







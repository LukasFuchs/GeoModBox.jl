# Momentum Equation

On geological time scales, Earth’s mantle and lithosphere do behave like a fluid and move and deform. A fluid does generally move due to forces acting on it whereas the forces must be in balance. In general, there are three major forces one might consider, i.e., inertia, surface, and volumetric forces. A common equation to describes such motion is given by: 

$\begin{equation}
\rho \frac{D \overrightharpoon{v}}{Dt} = \overrightharpoon{\nabla} \cdot \boldsymbol{\sigma} + \rho \boldsymbol{g},
\end{equation}$

where $\rho$ is the density [ $kg/m^3$ ], $\overrightharpoon{v}$ is the velocity vector [ $m/s$ ], $\boldsymbol{\sigma}$ is the *Cauchy stress tensor* [ $Pa$ ], $\boldsymbol{g}$ is the gravitational acceleration [ $m/s^2$ ], and the term on the left-hand side is the Lagrangian time derivative which is in Eulerian form $\frac{D}{Dt} = \frac{\partial{}}{\partial{t}} + \overrightharpoon{v} \cdot \overrightharpoon{\nabla}$. 

The *Cauchy stress tensor* is given by: 

$\begin{equation}
\boldsymbol{\sigma} = -\overrightharpoon{\nabla}{P} + \overrightharpoon{\nabla} \cdot \boldsymbol{\tau},
\end{equation}$

where $P$ is the total pressure ($P = P_{dynamic} + P_{hydrostatic}$) and $\boldsymbol{\tau}$ the *deviatoric stress tensor*. 

In Eulerian form, equation $(1)$ is given by (**Navier-Stokes equation**):

$\begin{equation}
\rho \left(\frac{\partial{v_{i}}}{\partial{t}} + v_{j}\frac{v_{i}}{\partial{x_{j}}}\right) = -\frac{\partial{P}}{\partial{x_{i}}} + \frac{\tau_{ij}}{\partial{x_j}} + \rho g_{i},
\end{equation}$

where summation over repeated indices is implied.

# Constitutive Relation

To solve equation $(3)$, one needs to define a rheology which, for a purely viscous medium, can be given by a constitutive relationship between stress and strain rate in the form of, e.g.:

$\begin{equation}
\tau_{ij} = 2 \eta \cdot \dot{\varepsilon}_{ij},
\end{equation}$

where $\eta$ is the dynamic viscosity in [ $Pa s$ ] and $\dot{\varepsilon}_{ij}$ the *strain rate tensor* in [ $1/s$ ] and given by: 

$\begin{equation}
\dot{\varepsilon}_{ij} = \frac{1}{2} \left(\frac{\partial{v_i}}{\partial{x_j}} + \frac{\partial{v_j}}{\partial{x_i}}\right),
\end{equation}$

# Stokes Equation

Assuming that the inertia forces are negligible in comparison to the gravitational forces, one can further simplify equation $(3)$ to:

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

To numerically solve for the three unknows $v_x$, $v_y$, and $P$, one needs to discretize the $x-$ and $y$-*component* of the **momentum equation** and the **continuum equation**. Here, we assume an incompressible medium, i.e. we use the so called Boussinesq-approximation. 

# Continuum Equation

Equation $(8)$ provides us two equations for our three unknowns. Thus, one needs to also consider the mass conservation equation (we do work with a continuum), where one can further simplify the problem by assuming an incompressible fluid (i.e., Boussinesq-approximation):

$\begin{equation}
\frac{\partial{v_i}}{\partial{x_i}} = 0.
\end{equation}$

Equations $(8)$ and $(9)$ enable us to solve for the three unknowns $v_x$, $v_y$, and $P$. 

# Equation of State

The buoyance term on the right-hand side of equation $(7)$, that is the density term which is temperature dependent (and pressure, but I do neglect this effect here so far), can be approximated with the so-called *equation of state* for the density. Here, its is a linear approximation of the change of density due to temperature variations and can be defined as:

$\begin{equation}
\rho = \rho_0 (1-\alpha T),
\end{equation}$

where $ρ_0$ is the reference density and $\alpha$ the thermal expansion coefficient [ $1/K$ ]. 

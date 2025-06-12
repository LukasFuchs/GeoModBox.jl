# Stokes Equation (1D)

Before tackling the *Stokes equation* in two dimensions, we start with a simpler, one-dimensional problem: uniaxial Stokes flow in a horizontal channel, assuming a known horizontal pressure gradient. This setup provides a first-order approximation for flows in magma or subduction channels. The one-dimensional Stokes equation in the $x$-direction is expressed as:

*$x$-component*

$\begin{equation}
0 = -\frac{\partial{P}}{\partial{x}} + \frac{\partial{\tau_{xy}}}{\partial{y}},
\end{equation}$

where $P$ is the pressure [Pa], $\frac{\partial}{\partial x_i}$ denotes the partial derivative in the $i$-th direction, and $\tau_{xy}$ is the horizontal shear stress [Pa], defined by:

$\begin{equation}
\tau_{xy} = 2 \eta \dot{\varepsilon}_{xy},
\end{equation}$

with $\eta$ the viscosity [Pa·s] and $\dot{\varepsilon}_{xy}$ the shear strain-rate [1/s], given by:

$\begin{equation}
\dot{\varepsilon}_{xy} = \frac{1}{2} \frac{\partial{v_x}}{\partial{y}}.
\end{equation}$

>**Note:** For the $y$-component of the Stokes equation, gravitational acceleration $g_y$ must be included.

## Discretization 

![Stokes1D_Grid](../assets/Stokes_1D_Grid.png)

**Figure 1.** **Channel flow setup and finite difference grid.** *Left:* Sketch of uniaxial channel flow driven by either a constant velocity at the top ($v_x$) and/or a horizontal pressure gradient $\left(\frac{\Delta P}{\Delta x} = P_1 - P_0\right)$, representing Couette, Poiseuille, or Couette-Poiseuille flow. *Right:* Finite difference grid with conservative gridding—viscosity is defined at *vertices*, and horizontal velocity is defined between them. The open circles at the top represent *ghost nodes* for horizontal velocity.

The finite difference grid employs a conservative scheme where velocity and viscosity are defined at staggered nodes. This ensures the horizontal shear stress is conserved across adjacent grid points, with stress values defined at the *vertices*. Such conservative gridding is essential when viscosity varies with depth. To build intuition, we first consider the simpler case of constant viscosity before addressing depth-dependent viscosity.

### Constant Viscosity

In the case of constant viscosity, a conservative gridding is not required. Equation (1) simplifies to:

$\begin{equation}
0 = -\frac{\partial{P}}{\partial{x}} + \eta\frac{\partial^2{v_x}}{\partial{y^2}}.
\end{equation}$

Applying a central difference approximation to the second derivative, and assuming a constant and known horizontal pressure gradient, this becomes:

$\begin{equation}
\frac{\partial{P}}{\partial{x}} = \eta \left( \frac{v_{x,j-1} - 2v_{x,j} + v_{x,j+1}}{\Delta{y^2}} \right),
\end{equation}$

which can be rewritten in a more compact form:

$\begin{equation}
\frac{\partial{P}}{\partial{x}}=av_{x,j-1}+bv_{x,j}+cv_{x,j+1}, 
\end{equation}$

where

$a = c = \frac{\eta}{\Delta{y^2}},\quad \textrm{and} \quad b = -\frac{2\eta}{\Delta{y^2}}.$

This results in a linear system of equations of the form:

$\begin{equation}
\bold{K} \cdot \vec{v_x} = \vec{rhs}
\end{equation}$ 

where $\bold{K}$ is a tridiagonal matrix, $\vec{v_x}$ is the vector of unknown horizontal velocites between the vertices, and $\vec{\text{rhs}}$ is the known right-hand side defined by the pressure gradient and boundary velocities. 

For simplicity, no separate solver for the constant viscosity case is included in `GeoModBox.jl`. Instead, viscosity is always treated numerically as an array, even in isoviscous cases. For implementation details, refer to the [source code](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/MomentumEquation/1Dsolvers.jl).

### Variable Viscosity

In the case of variable viscosity, Equation (1) becomes:

$\begin{equation}
0 = -\frac{\partial{P}}{\partial{x}} + \frac{\partial{\tau_{xy}}}{\partial{y}} = -\frac{\partial{P}}{\partial{x}} + \frac{\partial}{\partial{y}}\left(\eta\frac{\partial{v_x}}{\partial{y}}\right).
\end{equation}$

The differential operators in Equation (8) are approximated using central finite differences. The shear stress $\tau_{xy}$ and viscosity $\eta$ are defined at the vertices, while the velocity $v_x$ is defined at the *centroids* (the midpoints between adjacent vertices in 1D).

Using a central difference approximation for the shear stress, Equation (8) becomes:

$\begin{equation}
\frac{\partial{P}}{\partial{x}}=\frac{\tau_{xy,j+1}-\tau_{xy,j}}{\Delta{y}},\ \textrm{for}\ j = 1:nc, 
\end{equation}$

or, in terms of the velocity:

$\begin{equation}
\frac{\partial{P}}{\partial{x}}=\frac{\eta_{j+1}\frac{\partial{v_x}}{\partial{y}}\vert_{j+1}-\eta_{j}\frac{\partial{v_x}}{\partial{y}}\vert_{j}}{\Delta{y}},\ \textrm{for}\ j = 1:nc.
\end{equation}$

Approximating the velocity derivatives using central differences gives:

$\begin{equation}
\frac{\partial{P}}{\partial{x}}=\frac{\eta_{j+1}\frac{v_{x,j+1}-v_{x,j}}{\Delta{y}}-\eta_{j}\frac{v_{x,j}-v_{x,j-1}}{\Delta{y}}\vert_{j}}{\Delta{y}}.
\end{equation}$

>**Note**: The index $j$ ranges from $1$ to $nc$, where viscosity is defined on the vertices and velocity on the centroids.

Rewriting the above in terms of the unknown velocities, the equation becomes:

$\begin{equation}
\frac{\partial{P}}{\partial{x}}=av_{x,j-1}+bv_{x,j}+cv_{x,j+1}, 
\end{equation}$

where

$\begin{equation}
a = \frac{\eta_j}{\Delta{y^2}}, \quad 
b = -\frac{\eta_j+\eta_{j+1}}{\Delta{y^2}},\ \textrm{and}\ \quad 
c = \frac{\eta_{j+1}}{\Delta{y^2}}. 
\end{equation}$

This again results in a linear system of equations with a tridiagonal coefficient matrix.

### Boundary Conditions

To solve the equations, boundary conditions must be specified. For both *Dirichlet* and *Neumann* boundary conditions, the velocity values at the ghost nodes need to be defined. As with thermal boundary conditions, the ghost node velocities can be determined by assuming either a constant velocity at the boundary (Dirichlet) or a constant velocity gradient across the boundary (Neumann). The ghost node velocities are defined as:

**Dirichlet Boundary Conditions**

**Bottom**

$\begin{equation}
V_{G,S} = 2V_{BC,S} - v_{x,1}
\end{equation}$

**Top**

$\begin{equation}
V_{G,N} = 2V_{BC,N} - v_{x,nc}
\end{equation}$

**Neumann Boundary Conditions**

**Bottom**

$\begin{equation}
V_{G,S} = v_{x,1} - c_s\Delta{y},
\end{equation}$

**Top**

$\begin{equation}
V_{G,N}=v_{x,nc} + c_N\Delta{y},
\end{equation}$

where 

$\begin{equation}
c_S = \frac{dv_x}{dy}\vert_{S},\ \textrm{and}\ c_N=\frac{dv_x}{dy}\vert_{N}, 
\end{equation}$

are the velocity gradients across the boundary. 

To maintain a symmetric coefficient matrix, the coefficients at the centroids adjacent to the boundaries and the corresponding right-hand side (RHS) must be adjusted. The discretized equations at the bottom and top boundaries become:

**Dirichlet Boundary Conditions**

**Bottom**

$\begin{equation}
\left(b-a\right)v_{x,1}+cv_{x,2} = \frac{\partial{P}}{\partial{x}} - 2aV_{BC,S}
\end{equation}$

**Top**

$\begin{equation}
av_{x,nc-1}+\left(b-c\right)v_{x,nx} = \frac{\partial{P}}{\partial{x}} - 2cV_{BC,N}
\end{equation}$

**Neumann Boundary Conditions**

**Bottom**

$\begin{equation}
\left(b+a\right)v_{x,1}+cv_{x,2} = \frac{\partial{P}}{\partial{x}} + ac_SΔy
\end{equation}$

**Top**

$\begin{equation}
av_{x,nc-1}+\left(b+c\right)v_{x,nx} = \frac{∂P}{∂x} - cc_NΔy
\end{equation}$

### Solution 

There are different ways to solve the linear system of equations. The most convenient one is a direct solution using right division of the coefficient matrix by the right-hand side.

**Direct**

$\begin{equation}
v_x = \bold{K} ∖ rhs
\end{equation}$

Similar to the thermal problem, the system can also be solved using the **defect correction** method. This becomes especially useful when the system is non-linear, as it allows one to iteratively reduce the residual.

**Defect Correction**

The first step is to calculate the residual of the governing equation:

$\begin{equation}
R = -\frac{∂P}{∂x} + \frac{∂τ_{xy}}{∂y},
\end{equation}$

or in terms of the unknown horizontal velocity:

$\begin{equation}
R = -\frac{∂P}{∂x} + \bold{K} \cdot v_x. 
\end{equation}$

Assuming an initial guess for the horizontal velocity $v_{x,i}$, the initial residual is:

$\begin{equation}
R_i = -\frac{∂P}{∂x} + \bold{K_i} \cdot v_{x,i}.
\end{equation}$

Now assume that the initial guess can be corrected to the exact solution by adding a correction term $\delta{v_x}$. With some algebra:

$\begin{equation}
0 = -\frac{∂P}{∂x} + \bold{K}\left(v_{x,i}+ \delta{v_x} \right) = \bold{K}\cdot v_{x,i} -\frac{∂P}{∂x} + \bold{K}\cdot \delta{v_x} = R_i + \bold{K} \cdot \delta{v_x}.
\end{equation}$

Rearanging gives: 

$\begin{equation}
R_i = -\bold{K}\cdot{\delta{v_x}}, 
\end{equation}$

so the correction term is:

$\begin{equation}
\delta{v_x} = -\bold{K}^{-1}R_i. 
\end{equation}$

Finally, the corrected solution is:

$\begin{equation}
v_x^n = v_{x,i} + \delta{v_x}.
\end{equation}$

If the system is linear, this gives the solution in one step. For non-linear problems, this process must be iterated until the residual becomes sufficiently small.

For implementation details, see the [source code](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/MomentumEquation/1Dsolvers.jl).

An example of solving the channel flow using the defect correction method is provided in the [examples](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/StokesEquation/1D/ChannelFlow_1D.jl).

Solving the same problem using the direct method is part of the [exercises](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/08_1D_Stokes.ipynb).

# General Information

To numerically solve for the three unknows $v_x$, $v_y$, and $P$, one needs to discretize the x- and y-component of the **momentum equation** and the **continuum equation**. Here, we assume an incompressible medium, i.e. we use the so called Bousinesque approximation. 

---------------------
---------------------
## Stokes Equation (2D)

The Stokes equation (the **conservation equation of momentum**) in two dimensions is defined as: 

$$
\begin{equation}
0 = -\frac{\partial{P}}{\partial{x_i}} + \frac{\partial}{\partial{x_j}}\tau_{ij} + \rho \cdot g_i, 
\end{equation}
$$

where $P$ is the total pressure,  $\rho$ is the density, $g_i$ the gravitational acceleration, $\frac{\partial}{\partial_i}$ the spatial derivative in the direction of $x_i$, and $\tau_{ij}$ is the deviatoric stress tensor and defined as: 

$$\begin{equation}
\tau_{ij} = 2\eta \dot{\varepsilon}_{ij}.
\end{equation}
$$

The **conservation equation of mass** is defined as: 

$$\begin{equation}
div\left(\overrightarrow{v} \right) = \left(\frac{\partial{v_i}}{\partial{x_i}}+\frac{\partial{v_j}}{\partial{x_j}}\right) = 0, 
\end{equation}
$$

where $v_i$ is the velocity in the i-th direction. 

To approximate the partial derivative operators using the finite difference formulations, one needs to first define the numerical grid of the domain and the location of the corresponding parameters.

### Discretization 

The conservation equations of *momentum* and *mass* are solved properly in two dimensions (*x* and *y*) using a staggered finite difference grid, where the horizontal (*cyan dashes*) and vertical (*orange dashes*) velocity are defined in between the regular grid points or *vertices*, and the pressure (*red circles*) within a finite difference cells or *centroids* (Figure 1). A staggered grid enables the conservation of the stress between adjacent grid points and one can solve equations $(1)$ and $(3)$ for the unknows.  

<img src="./Figures/MomentumGrid_2D.png" alt="drawing" width="600"/> <br>
**Figure 1. Staggered finite difference grid.** Discretization of the conservation equations of momemtum and mass. The horizontal and vertical velocities require *ghost nodes* at the north and south and east and west boundary, respectively. <!-- What about the pressure ghost nodes? Are they even necessary? -->

### Constant Viscosity

Let's first assume a special case of the Stokes equation, a **constant viscosity**, which simplifies equation $(1)$ to: 

$$\begin{equation}
0 = -\frac{\partial{P}}{\partial{x_i}} + 2\eta\frac{\partial^2{v_i}}{\partial{x_i^2}} + \eta\left(\frac{\partial^2{v_i}}{\partial{x_j^2}}+\frac{\partial^2{v_j}}{\partial{x_i^2}}\right) + \rho \cdot g_i. 
\end{equation}
$$

Using equation $(3)$ and neglecting a horizontal graviational acceleration, equation $(4)$ can further be simplified to: 

*x-component* 

$$
\begin{equation}
-\frac{\partial{P}}{\partial{x}} + \eta\frac{\partial^2{v_x}}{\partial{x^2}} + \eta\frac{\partial^2{v_x}}{\partial{y^2}} = 0, 
\end{equation}
$$

*y-component*

$$
\begin{equation}
-\frac{\partial{P}}{\partial{y}} + \eta\frac{\partial^2{v_y}}{\partial{y^2}} + \eta\frac{\partial^2{v_y}}{\partial{x^2}} = - \rho g_y. 
\end{equation}
$$

For each equation, one can define a so-called numerical stencil which highlights the position of the parameters with respect to a central point $(i,j)$, where *i* and *j* are the indeces in the horizontal and vertical direction, respectively. The central point corresponds to the number of the equation in the linear system of equations. 

#### Stencil

The stencils for the momentum equation assuming a constant viscosity shows the points required to solve the equation using the finite difference approach (Figure 2). 

<img src="./Figures/Stencil_const_eta.png" alt="drawing" width="600"/> <br>
**Figure 2. Numerical stencils for the momentum equation.** a) *x-component. b) y-component.

Using the finite difference operators equations $(5)$ and $(6)$ are defined as: 

*x-component*

$$\begin{equation}
-\frac{P_{i,j}-P_{i-1,j}}{\Delta{x}} + \eta\frac{v_{x,(i-1,j)}-2v_{x,(i,j)}+v_{x,(i+1,j)}}{\Delta{x^2}} + \eta\frac{v_{x,(i,j-1)}-2v_{x,(i,j)}+v_{x,(i,j+1)}}{\Delta{y^2}} = 0
\end{equation}
$$

$$
\begin{equation}
P_CP_{i,j} + P_WP_{i-1,j} + Sv_{x,(i,j-1)} +  Wv_{x,(i-1,j)} + Cv_{x,(i,j)} + E v_{x,(i+1,j)} + N v_{x,(i,j+1)} = 0, 
\end{equation}
$$

where 

$$\begin{equation}
\begin{split}
P_C = -\frac{1}{\Delta{x}},  \\
P_W = \frac{1}{\Delta{x}}, \\
S = \frac{\eta}{\Delta{y^2}}, \\
W = \frac{\eta}{\Delta{x^2}}, \\
C = -2\eta\left(\frac{1}{\Delta{x^2}}+\frac{1}{\Delta{y^2}}\right), \\
E = \frac{\eta}{\Delta{x^2}}, \\
N = \frac{\eta}{\Delta{y^2}}. \\
\end{split}
\end{equation}
$$

*y-component*

$$\begin{equation}
-\frac{P_{i,j}-P_{i,j-1}}{\Delta{y}} + \eta\frac{v_{y,(i,j-1)}-2v_{y,(i,j)}+v_{x,(i,j+1)}}{\Delta{y^2}} + \eta\frac{v_{y,(i-1,j)}-2v_{y,(i,j)}+v_{y,(i+1,j)}}{\Delta{x^2}} = -\rho g_y
\end{equation}
$$

$$
\begin{equation}
P_SP_{i,j-1} + P_CP_{i,j} + Sv_{y,(i,j-1)} + Wv_{y,(i-1,j)} + Cv_{y,(i,j)} + E v_{y,(i+1,j)} + N v_{y,(i,j+1)} = -\rho g_y, 
\end{equation}
$$

where 

$$\begin{equation}
\begin{split}
P_S = \frac{1}{\Delta{y}},  \\
P_C = -\frac{1}{\Delta{y}}, \\
S = \frac{\eta}{\Delta{y^2}}, \\
W = \frac{\eta}{\Delta{x^2}}, \\
C = -2\eta\left(\frac{1}{\Delta{x^2}}+\frac{1}{\Delta{y^2}}\right), \\
E = \frac{\eta}{\Delta{x^2}}, \\
N = \frac{\eta}{\Delta{y^2}}. \\
\end{split}
\end{equation}
$$

#### Boundary Conditions

The most common boundary conditions for the momentum equation are a combination of *Dirichlet* and *Neumann* velocity boundary conditions: **free slip** and **no slip** boundary condition. 

---------------------

**Free slip** 

Free slip boundary conditions allow the fluid to move along the boundary assuming **no shear stress** and **no orthogonal velocity**. That is for the lateral boundaries **(E, W)** the conditions are: 

$$\begin{equation}
\begin{split}
v_x = 0, \\
\frac{\partial{v_{y}}}{\partial{x}}=0,
\end{split}
\end{equation}
$$

and for the horizontal bounaries **(N, S)** the conditions are: 

$$\begin{equation}
\begin{split}
v_y = 0, \\
\frac{\partial{v_{x}}}{\partial{y}}=0.
\end{split}
\end{equation}
$$

Solving the system of equations using a direct solution, one needs to modify the coefficients of the required points adjacent to the boundaries using the velocities defined at the *ghost nodes* **and** the right-hand side. 

*x-component*

For free slip boundary conditions one needs to define the horizontal velocity $v_x$ for the *ghost nodes* at the upper (**N**) and lower (**S**) boundary which are defined as: 

*South* (j = 1)

$$\begin{equation}
v_{x,(i,GS)} = v_{x,(i,1)},
\end{equation}
$$

*North* (j = ncy)

$$\begin{equation}
v_{x,(i,GN)} = v_{x,(i,ncy)}.
\end{equation}
$$

Along the lateral boundaries (**E**, **W**) we can simply set the velocities to zero.

*y-component*

For the vertical velocities $v_y$ one needs to define the velocity for the *ghost nodes* at the left (**W**) and right (**E**) boundary as: 

*West* (i = 1)

$$\begin{equation}
v_{y,(GW,j)} = v_{y,(1,j)},
\end{equation}
$$

*East* (i = ncx)

$$\begin{equation}
v_{y,(GE,j)} = v_{y,(ncx,j)}.
\end{equation}
$$

Along the horizontal boundaries (**N**, **S**), we can simply set the velocity $v_y$ to zero. 

Using the equations $(15)$ - $(18)$ the coefficients of the equations adjacent to the corresponding boundaries changes to (the right-hand side actually does not need to be changed for free slip boundary conditions): 

*x-component* 

*South* (j = 1)

$$\begin{equation}
P_WP_{i-1,j}+P_CP_{i,j}+Wv_{x,(i-1,j)}+Cv_{x,(i,j)}+Ev_{x,(i+1,j)}+Nv_{x,(i,j+1)} = 0,
\end{equation}
$$

where $C = -\frac{2\eta}{\Delta{x^2}}-\frac{\eta}{\Delta{y^2}}$.

*North* (j = ncy)

$$\begin{equation}
P_WP_{i-1,j}+P_CP_{i,j}+Sv_{x,(i,j-1)}+Wv_{x,(i-1,j)}+Cv_{x,(i,j)}+Ev_{x,(i+1,j)} = 0,
\end{equation}
$$

where $C = -\frac{2\eta}{\Delta{x^2}}-\frac{\eta}{\Delta{y^2}}$.

Along the lateral boundaries $C=1$ and the remaining coefficients are equal to zero.

*y-component* 

*West* (i = 1)

$$\begin{equation}
P_SP_{i,j-1}+P_CP_{i,j}+Sv_{y,(i,j-1)}+Cv_{y,(i,j)}+Ev_{y,(i+1,j)}+Nv_{y,(i,j+1)} = -\rho g_y,
\end{equation}
$$

where $C = -\frac{\eta}{\Delta{x^2}}-\frac{2\eta}{\Delta{y^2}}$.

*East* (i = ncx)

$$\begin{equation}
P_SP_{i,j-1}+P_CP_{i,j}+Sv_{y,(i,j-1)}+Wv_{y,(i-1,j)}+Cv_{y,(i,j)}+Nv_{y,(i,j+1)} = -\rho g_y,
\end{equation}
$$

where $C = -\frac{\eta}{\Delta{x^2}}-\frac{2\eta}{\Delta{y^2}}$.

Along the horizontal boundaries $C=1$ and the remaining coefficients are equal to zero.

--------------

**No slip** 

No slip boundary conditions fix the fluid along the boundary and sets the horizontal ($v_x$) and vertical ($v_y$) velocity equal to zero.. 

That is, for all boundaries **(E, W, S, N)** the conditions are: 

$$\begin{equation}
\begin{split}
v_x = 0, \\
v_y = 0.
\end{split}
\end{equation}
$$

*x-component*

One needs the velocity at the *ghost nodes* for the horizontal velocity $v_x$ at the lower (**S**) and upper (**N**) boundary which are defined as: 

*South* (j = 1)

$$\begin{equation}
v_{x,(i,GS)} = 2V_{BC,S} - v_{x,(i,1)},
\end{equation}
$$

where $V_{BC,S}$ is the velocity along the boundary (here 0).

*North* (j = ncy)

$$\begin{equation}
v_{x,(i,GN)} = 2V_{BC,N} - v_{x,(i,ncy)},
\end{equation}
$$

where $V_{BC,N}$ is the velocity along the boundary (here 0). Along the lateral boundaries (**E**, **W**) we can simply set the velocity $v_x$ to zero. 

*y-component* 

For the vertical velocities $v_y$ one needs to define the velocity for the *ghost nodes* at the left (**W**) and right (**E**) boundary as: 

*West* (i = 1)

$$\begin{equation}
v_{y,(GW,j)} = 2V_{BC,W} - v_{y,(1,j)},
\end{equation}
$$

where $V_{BC,W}$ is the velocity along the boundary (here 0).

*East* (i = ncx)

$$\begin{equation}
v_{y,(GE,j)} = 2V_{BC,E} - v_{y,(ncx,j)},
\end{equation}
$$

where $V_{BC,E}$ is the velocity along the boundary (here 0). Along the horizontal boundaries (**N**, **S**), we can simply set the velocity $v_y$ to zero. 

Using the equations $(24)$ - $(27)$ the coefficients of the equations adjacent to the corresponding boundaries and the right-hand side changes to: 

*x-component* 

*South* (j = 1)

$$\begin{equation}
P_WP_{i-1,j}+P_CP_{i,j}+Wv_{x,(i-1,j)}+Cv_{x,(i,j)}+Ev_{x,(i+1,j)}+Nv_{x,(i,j+1)} = -2\frac{\eta}{\Delta{y^2}}V_{BC,S},
\end{equation}
$$

where $C = -\frac{2\eta}{\Delta{x^2}}-\frac{3\eta}{\Delta{y^2}}$.

*North* (j = ncy)

$$\begin{equation}
P_WP_{i-1,j}+P_CP_{i,j}+Sv_{x,(i,j-1)}+Wv_{x,(i-1,j)}+Cv_{x,(i,j)}+Ev_{x,(i+1,j)} = -2\frac{\eta}{\Delta{y^2}}V_{BC,N},
\end{equation}
$$

where $C = -\frac{2\eta}{\Delta{x^2}}-\frac{3\eta}{\Delta{y^2}}$.

Along the lateral boundaries $C=1$ and the remaining coefficients and the right-hand side are equal to zero.

*y-component* 

*West* (i = 1)

$$\begin{equation}
P_SP_{i,j-1}+P_CP_{i,j}+Sv_{y,(i,j-1)}+Cv_{y,(i,j)}+Ev_{y,(i+1,j)}+Nv_{y,(i,j+1)} = -\rho g_y - 2\frac{\eta}{\Delta{x^2}}V_{BC,W},
\end{equation}
$$

where $C = -\frac{3\eta}{\Delta{x^2}}-\frac{2\eta}{\Delta{y^2}}$.

*East* (i = ncx)

$$\begin{equation}
P_SP_{i,j-1}+P_CP_{i,j}+Sv_{y,(i,j-1)}+Wv_{y,(i-1,j)}+Cv_{y,(i,j)}+Nv_{y,(i,j+1)} = -\rho g_y - 2\frac{\eta}{\Delta{x^2}}V_{BC,W},
\end{equation}
$$

where $C = -\frac{3\eta}{\Delta{x^2}}-\frac{2\eta}{\Delta{y^2}}$.

Along the horizontal boundaries $C=1$ and the remaining coefficients and the right-hand side are equal to zero.

### Variable Viscosity

#### Stencil

#### Boundary Conditions

### Continuum Equation 

#### Stencil

## Solution 

### Direct 

### Defect Correction

-----------------------
-----------------------

## Examples

<!-- 
- 2D case 
-- Discretized equations
-- Solving the equations 
--- Direct solution 
--- Defection corrections solution
-->

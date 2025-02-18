# General Information (Energyequation 2D)

&emsp; This directory provides several examples solving the diffusive part of the *temperature conservation equation* (in 2-D, steady state and time-dependent) using different numerical discretization methods. The different solvers for a 2-D problem are located in [src/HeatEquation/](../../../src/HeatEquation/).

## Discretization 

...

<img src="./Figures/Diffusion_Grid.jpg" alt="drawing" width="600"/> <br>
**Figure 1.** ... 

&emsp; In two dimesions (*x* and *y*), the diffusive part of *temperature conservation equation* is described by (assuming only radiogenic heating):

$$
\begin{equation}
\rho c_p \frac{\partial T}{\partial t} = -\frac{\partial q_x}{\partial x} -\frac{\partial q_y}{\partial y} + \rho H, 
\end{equation}
$$

or including Fourier’s law (assuming variable thermal parameters):

$$
\begin{equation}
\rho c_p \frac{\partial T}{\partial t} = \frac{\partial}{\partial x} k \frac{\partial T}{\partial x} + \frac{\partial}{\partial y} k \frac{\partial T}{\partial y} + \rho H.
\end{equation}
$$

Assuming that the thermal parameters are constant, equation $(2)$ simplifies to: 

$$
\begin{equation}
\frac{\partial T}{\partial t} = \kappa (\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2}) + \frac{Q}{\rho c_p},
\end{equation}
$$
  
where $\kappa = k/\rho /c_p$ is the thermal diffusivity [m<sup>2</sup>/s] and $Q=\rho H$ is the heat production rate per volume [W/m<sup>3</sup>]. Equation $(3)$ is a *parabolic partial differential equation* which can be solved numerically in different manners, assuming initial and boundary conditions are defined. 

&emsp;First we would like to discuss a simple, but effective, finite difference method to discretize and solve the equation, that is the forward in time and centered in space (FTCS) method in an *explicit* manner. This finite difference scheme will converge to the exact solution for small $\Delta x$, $\Delta y$, and $\Delta t$. The advantage of an explicit description is that it is **simple** to derive and rather **fast** computationally. However, it is only numerically stable as long as the *heat diffusive stability criterion* is fulfilled. The stability criterion can be determined by a *Von Neumann* stability analysis, which analyzes the growth of an eigenmode perturbation for a certain finite difference approach. In case of an **explicit 2-D finite difference approach**, the heat diffusive stability criterion is defined as 
$$ \begin{equation}
\Delta t < \frac{minimum({\Delta x^2,\Delta y^2})}{3 \kappa}
\end{equation}
$$
(assuming equal grid spacing in *x*- and *y*-direction). Attention, the stability criterion limits the model time step by the model’s resolution. 

### Explicit, FTCS (or Forward Euler Method)

...

$$
\begin{equation}
\frac{T_{i,j}^{n+1} - T_{i,j}^{n} }{\Delta t} = \kappa \left( \frac{T_{i-1,j}^{n} - 2T_{i,j}^{n} + T_{i+1,j}^{n}}{\Delta{x}^2} + \frac{T_{i,j-1}^{n} - 2T_{i,j}^{n} + T_{i,j+1}^{n}}{\Delta{z^2}} \right) + \frac{Q_{i,j}^n}{\rho c_p}
\end{equation}
$$

$$
\begin{equation}
T_{i,j}^{n+1} = T_{i,j}^{n} + a\left(T_{i-1,j}^{n} - 2T_{i,j}^{n} + T_{i+1,j}^{n}\right) + b\left(T_{i,j-1}^{n} - 2T_{i,j}^{n} + T_{i,j+1}^{n}\right) + \frac{Q_{i,j}^n \Delta{t}}{\rho c_p}, 
\end{equation}
$$

where $a = \frac{\kappa \Delta{t}}{\Delta{x^2}}$ and $b = \frac{\kappa \Delta{t}}{\Delta{y^2}}$. 

...


#### Boundary Conditions

&emsp;Different thermal boundary conditions can be set for our model for which we now utilize the *ghost nodes*. Here, we focus on two fundamental conditions, the *Dirichlet* and *Neumann* boundary conditions. To consider both boundary conditions to solve the equations, one needs to define the temperature at the *ghost nodes* based on the assumed boundary condition. The Dirichlet boundary condition defines a constant temperature along the boundary, such that the temperatures at the left (west), right (east), bottom (south), and top (north) *ghost nodes* are defined as:

$$
\begin{equation}
T_{G,W} = 2T_{BC,W} - T_{1,:},
\end{equation}
$$

$$
\begin{equation}
T_{G,E} = 2T_{BC,E} - T_{ncx,:},
\end{equation}
$$

$$
\begin{equation}
T_{G,S} = 2T_{BC,S} - T_{:,1},
\end{equation}
$$

$$
\begin{equation}
T_{G,N} = 2T_{BC,N} - T_{:,ncy},
\end{equation}
$$

where $T_{G,W}$, $T_{G,E}$, $T_{G,S}$, and $T_{G,N}$ and $T_{BC,W}$, $T_{BC,E}$, $T_{BC,S}$, and $T_{BC,N}$ are the temperature at the left, right, bottom, and top *ghost nodes* and the constant temperatures at the left, right, bottom, and top boundary, respectrively. 

&emsp;The Neumann boundary condition defines that the variation of a certain parameter across the boundary is defined, that is, for example, the temperature across the boundary or thermal heat flux *q* through the boundary. The temperature at the *ghost nodes* is then defined as: 

$$
\begin{equation}
T_{G,W} = T_{1,:} - c_{W} \Delta{x},
\end{equation}
$$

$$
\begin{equation}
T_{G,E} = T_{ncx,:} + c_{E} \Delta{x},
\end{equation}
$$

$$
\begin{equation}
T_{G,S} = T_{:,1} - c_{S} \Delta{y},
\end{equation}
$$

$$
\begin{equation}
T_{G,N} = T_{:,ncy} + c_{N} \Delta{y},
\end{equation}
$$

where 

$$
\begin{equation}
\left. c_{W} = \frac{\partial{T}}{\partial{x}} \right\vert_{W}, \left. c_{E} = \frac{\partial{T}}{\partial{x}} \right\vert_{E}, 
\left. c_{S} = \frac{\partial{T}}{\partial{y}} \right\vert_{S},
\left. c_{N} = \frac{\partial{T}}{\partial{y}} \right\vert_{N},
\end{equation}
$$

are the constant heat fluxes along the left, right, bottom, and top boundary, respectively. 

&emsp; Now one can solve equation $(6)$ for temperature on the centroids at the next time step using the defined temperature at the *ghost nodes* if necessary.  

## Numerical Schemes

...

### Implicit, FTCS (or Backward Euler Method)

&emsp;The fully implicit finite difference scheme is unconditionally stable and one can use time steps larger than the diffusion time criterion. In 2-D, the *temperature conservation equation* is then given as: 

$$
\begin{equation}
\frac{T_{i,j}^{n+1}-T_{i,j}^n}{\Delta t} = 
\kappa \left( 
    \frac{T_{i-1,j}^{n+1}-2T_{i,j}^{n+1}+T_{i+1,j}^{n+1}}{\Delta x^2} + 
    \frac{T_{i,j-1}^{n+1}-2T_{i,j}^{n+1}+T_{i,j+1}^{n+1}}{\Delta y^2} 
    \right) + 
\frac{Q_{i}^n}{\rho c_p},
\end{equation}
$$

where *n* is the current and *n+1* the next time step, $\Delta{t}$ is the time step length, $\Delta{x}$ and $\Delta{y}$ are the horizontal and vertical grid spacing, and *i* and *j* are the horizontal and vertical indeces, respectively. Rearranging equation $(15)$ into known and unknown variables, one obtains a linear system of equations in the form of: 

$$
\begin{equation}
-b T_{i,j-1}^{n+1} - a T_{i-1,j}^{n+1} + 
\left(2a + 2b + c \right) T_{i,j}^{n+1} - 
a T_{i+1,j}^{n+1} - b T_{i,j+1}^{n+1} = 
c T_{i,j}^n + \frac{Q_{i,j}^n}{\rho c_p},
\end{equation}
$$

where $a = \frac{\kappa}{\Delta{x^2}}$, $b = \frac{\kappa}{\Delta{y^2}}$, and $c = \frac{1}{\Delta{t}}$. This is a linear system of equation in the form of $\boldsymbol{A}\cdot x = rhs$, where $\boldsymbol{A}$ is a coefficient matrix (here with five non-zero diagonals), $x$ the unknown vector, and $rhs$ the known right-hand side. The main advantage of the implicit method is that there are no restrictions on the time step, but this does not mean that it is more accurate. Taking too large time steps may result in an inaccurate solution for features with small spatial scales.

#### Boundary Conditions

&emsp;The temperature on the *ghost nodes* to solve the equations on the centroids adjacent to the boundary are defined as before (equations $(6)$ and $(13)$). To obtain a symmetric coefficient matrix to solve the linear system of euqations, however, one needs to modify the coefficients of the centroids adjacent to the boundary and the corresponding right-hand side, such that the equations are defined as:  

**Dirichlet** <br>
*West (i=1)*

$$
\begin{equation}
-b T_{1,j-1}^{n+1} + 
\left(3 a + 2b + c\right) T_{1,j}^{n+1} 
- a T_{2,j}^{n+1}  - b T_{1,j+1}^{n+1} = 
c T_{1,j}^{n} + 2 a T_{BC,W} + \frac{Q_{i,j}}{\rho c_p}, 
\end{equation}
$$

*East (i = ncx)*

$$
\begin{equation}
- bT_{ncx,j-1}^{n+1} - aT_{ncx-1,j}^{n+1} + 
\left(3 a + 2b + c\right) T_{ncx,j}^{n+1} 
- b T_{ncx,j+1}^{n+1} = 
c T_{ncx,j}^{n} + 2 a T_{BC,E} + \frac{Q_{i,j}}{\rho c_p}, 
\end{equation}
$$

*South (j = 1)*

$$
\begin{equation}
- aT_{i-1,1}^{n+1} + 
\left(2a + 3b + c\right) T_{i,1}^{n+1} 
- a T_{i+1,1}^{n+1} - bT_{i,2}^{n+1} = 
c T_{i,1}^{n} + 2 b T_{BC,S} + \frac{Q_{i,j}}{\rho c_p}, 
\end{equation}
$$

*North (j = ncy)*

$$
\begin{equation}
- bT_{i,ncy}^{n+1} - aT_{i-1,ncy}^{n+1} + 
\left(2a + 3b + c\right) T_{i,ncy}^{n+1} 
- a T_{i+1,ncy}^{n+1} = 
c T_{i,ncy}^{n} + 2 b T_{BC,N} + \frac{Q_{i,j}}{\rho c_p}, 
\end{equation}
$$

**Neumann** <br>
*West (i=1)*

$$
\begin{equation}
-b T_{1,j-1}^{n+1} + 
\left(a + 2b + c\right) T_{1,j}^{n+1} 
- a T_{2,j}^{n+1}  - b T_{1,j+1}^{n+1} = 
c T_{1,j}^{n} - a c_W \Delta{x} + \frac{Q_{i,j}}{\rho c_p}, 
\end{equation}
$$

*East (i = ncx)*

$$
\begin{equation}
- bT_{ncx,j-1}^{n+1} - aT_{ncx-1,j}^{n+1} + 
\left(a + 2b + c\right) T_{ncx,j}^{n+1} 
- b T_{ncx,j+1}^{n+1} = 
c T_{ncx,j}^{n} + a c_E \Delta{x} + \frac{Q_{i,j}}{\rho c_p}, 
\end{equation}
$$

*South (j = 1)*

$$
\begin{equation}
- aT_{i-1,1}^{n+1} + 
\left(2a + b + c\right) T_{i,1}^{n+1} 
- a T_{i+1,1}^{n+1} - bT_{i,2}^{n+1} = 
c T_{i,1}^{n} - b c_S \Delta{y} + \frac{Q_{i,j}}{\rho c_p}, 
\end{equation}
$$

*North (j = ncy)*

$$
\begin{equation}
- bT_{i,ncy}^{n+1} - aT_{i-1,ncy}^{n+1} + 
\left(2a + b + c\right) T_{i,ncy}^{n+1} 
- a T_{i+1,ncy}^{n+1} = 
c T_{i,ncy}^{n} + 2 c_N \Delta{y} + \frac{Q_{i,j}}{\rho c_p}, 
\end{equation}
$$

### Defection Correction Method

&emsp; The defection correction method is an iterative solution, in which the residual of the diffusion equation for an initial temperature conditions is reduced by a correction term. In case, the system is linear, one iteration is sufficient enough to optain the exact solution. 

*Theory*

The diffusion equation, in an implicit form, can be simplified to an equation in the form of: 

$$
\begin{equation}
\boldsymbol{K} \cdot T - b = R, 
\end{equation}
$$

where $\boldsymbol{K}$ is the coefficient matrix, $T$ is the temperature at the new time step, $b$ is an term containing the remaining variables, and $R$ is the resiual. Assuming an initial temperature guess $T_i$, the initial residual $R_i$ is given by: 

$$
\begin{equation}
R_i = \boldsymbol{K} \cdot T_i - b.
\end{equation}
$$

Adding a correction term $\delta{T}$ to the initial guess, assuming that it results in zero residual, leads to: 

$$
\begin{equation}
0 = \boldsymbol{K} \left(T_i + \delta{T} \right) - b = \boldsymbol{K} T_i - b + \boldsymbol{K} \delta{T} = R_i + \boldsymbol{K} \delta{T},
\end{equation}
$$

which results in:

$$
\begin{equation}
R_i = -\boldsymbol{K} \delta{T}, 
\end{equation}
$$

and finally the correction term: 

$$
\begin{equation}
\delta{T} = -\boldsymbol{K}^{-1} R_i. 
\end{equation}
$$

The coefficients of the matrix can be derived, for example, via: 

$$
\begin{equation}
\frac{\partial{T}}{\partial{t}} - \kappa \left( \frac{\partial^2{T}}{\partial{x}^2} + \frac{\partial^2{T}}{\partial{y}^2} \right) - \frac{Q_{i,j}^n}{\rho c_p} = R, 
\end{equation}
$$

$$
\begin{equation}
\frac{T_i^{n+1}-T_i^{n}}{\Delta{t}} - \kappa 
\left( \frac{T_{i-1,j}^{n+1} - 2 T_{i,j}^{n+1} + T_{i+1,j}^{n+1}}{\Delta{x}^2} + \frac{T_{i,j-1}^{n+1} - 2 T_{i,j}^{n+1} + T_{i,j+1}^{n+1}}{\Delta{y}^2}  
\right) - \frac{Q_{i,j}^n}{\rho c_p} = R,
\end{equation}
$$

$$
\begin{equation}
-b T_{i,j-1}^{n+1} - a T_{i-1,j}^{n+1} + 
\left(2a + 2b + c \right) T_{i,j}^{n+1} - 
a T_{i+1,j}^{n+1} - b T_{i,j+1}^{n+1} - 
c T_{i,j}^n - \frac{Q_{i,j}^n}{\rho c_p} = 
R,
\end{equation}
$$

where

$$
\begin{equation}
a = \frac{\kappa}{\Delta{x}^2},\  
b = \frac{\kappa}{\Delta{y}^2},\ 
\ and \ c = \frac{1}{\Delta{t}}, 
\end{equation}
$$

&emsp;For more details on how this is implemented, see [*BackwardEuler.jl*](../../../src/HeatEquation/BackwardEuler.jl).

### Cranck-Nicolson Approach (CNA)

&emsp;The fully implicit method works well, but is only first order accurate in time. A way to modify this is to employ a Crank-Nicolson time step discretization, which is implicit and thus second order accurate in time. In 2-D, equation $(3)$ is then described as: 

$$
\begin{equation}
\frac{T_{i,j}^{n+1} - T_{i,j}^{n}}{\Delta t} = 
\frac{\kappa}{2}\frac{(T_{i-1,j}^{n+1}-2T_{i,j}^{n+1}+T_{i+1,j}^{n+1})+(T_{i-1,j}^{n}-2T_{i,j}^{n}+T_{i+1,j}^{n})}{\Delta x^2} + 
\frac{\kappa}{2}\frac{(T_{i,j-1}^{n+1}-2T_{i,j}^{n+1}+T_{i,j+1}^{n+1})+(T_{i,j-1}^{n}-2T_{i,j}^{n}+T_{i,j+1}^{n})}{\Delta y^2}
\end{equation}
$$

&emsp;Similar to the implicit method, we need to modify the coefficients and the right-hand side using different boundary conditions to obtain a symmetric coefficient matrix. Thus, the equations for the centroids adjacent to the boundaries are defined as: 

**Dirichlet**<br>
*West (i=1)*

$$
\begin{equation}
- b T_{1,j-1}^{n+1} +
\left(3a + 2b + c \right) T_{1,j}^{n+1} 
- a T_{2,j}^{n+1} - b T_{1,j+1}^{n+1}
= 
b T_{1,j-1}^{n} -
\left( 3a + 2b - c \right) T_{1,j}^{n} 
+ a T_{2,j}^{n} 
+ b T_{1,j+1}^{n}
+ 4 a T_{BC,W},
\end{equation}
$$

*East (i = ncx)*

$$
\begin{equation}
- b T_{ncx,j-1}^{n+1} - a T_{ncx-1,j}^{n+1} +
\left(3a + 2b + c \right) T_{ncx,j}^{n+1} 
- b T_{ncx,j+1}^{n+1}
= 
b T_{ncx,j-1}^{n} + a T_{ncx-1,j}^{n} -
\left( 3a + 2b - c \right) T_{ncx,j}^{n} 
+ b T_{ncx,j+1}^{n}
+ 4 a T_{BC,E},
\end{equation}
$$

*South (j = 1)*

$$
\begin{equation}
- a T_{i-1,1}^{n+1} +
\left(2a + 3b + c \right) T_{i,1}^{n+1} 
- a T_{i+1,1}^{n+1}
- b T_{i,2}^{n+1}
= 
a T_{i-1,1}^{n} -
\left( 2a + 3b - c \right) T_{i,1}^{n} 
+ a T_{i+1,1}^{n}
+ b T_{i,2}^{n}
+ 4 b T_{BC,S},
\end{equation}
$$

*North (j = ncy)*

$$
\begin{equation}
- b T_{i,ncy-1}^{n+1} +
- a T_{i-1,ncy}^{n+1} +
\left(2a + 3b + c \right) T_{i,ncy}^{n+1} 
- a T_{i+1,ncy}^{n+1}
= 
b T_{i,ncy-1}^{n} +
a T_{i-1,ncy}^{n} -
\left( 2a + 3b - c \right) T_{i,ncy}^{n} 
+ a T_{i+1,ncy}^{n}
+ 4 b T_{BC,N}.
\end{equation}
$$


**Neumann**<br>

*West (i=1)*

$$
\begin{equation}
- b T_{1,j-1}^{n+1} +
\left(a + 2b + c \right) T_{1,j}^{n+1} 
- a T_{2,j}^{n+1} - b T_{1,j+1}^{n+1}
= 
b T_{1,j-1}^{n} -
\left( a + 2b - c \right) T_{1,j}^{n} 
+ a T_{2,j}^{n} 
+ b T_{1,j+1}^{n}
- 2 a c_W \Delta{x},
\end{equation}
$$

*East (i = ncx)*

$$
\begin{equation}
- b T_{ncx,j-1}^{n+1} - a T_{ncx-1,j}^{n+1} +
\left(a + 2b + c \right) T_{ncx,j}^{n+1} 
- b T_{ncx,j+1}^{n+1}
= 
b T_{ncx,j-1}^{n} + a T_{ncx-1,j}^{n} -
\left( a + 2b - c \right) T_{ncx,j}^{n} 
+ b T_{ncx,j+1}^{n}
+ 2 a c_E \Delta{x},
\end{equation}
$$

*South (j = 1)*

$$
\begin{equation}
- a T_{i-1,1}^{n+1} +
\left(2a + b + c \right) T_{i,1}^{n+1} 
- a T_{i+1,1}^{n+1}
- b T_{i,2}^{n+1}
= 
a T_{i-1,1}^{n} -
\left( 2a + b - c \right) T_{i,1}^{n} 
+ a T_{i+1,1}^{n}
+ b T_{i,2}^{n}
- 2 b c_S \Delta{y}
\end{equation}
$$

*North (j = ncy)*

$$
\begin{equation}
- b T_{i,ncy-1}^{n+1} +
- a T_{i-1,ncy}^{n+1} +
\left(2a + b + c \right) T_{i,ncy}^{n+1} 
- a T_{i+1,ncy}^{n+1}
= 
b T_{i,ncy-1}^{n} +
a T_{i-1,ncy}^{n} -
\left( 2a + b - c \right) T_{i,ncy}^{n} 
+ a T_{i+1,ncy}^{n}
+ 2 b c_N \Delta{y}
\end{equation}
$$

&emsp;However, the band-width of the coefficient matrix increases as in the fully implicit case. Thus, the method becomes memory intensiv for models with a high resoltuion. For more details on how this is implemented, see [*CNA.jl*](../../../src/HeatEquation/CNA.jl).

### Alternating Direction Implicit
&emsp; Within the ADI method, one basically decomposes the calculation of one time step into two half-steps. For the first step $(n \ \text{to}\ n+1/2)$, the *x*-direction is solved explicitly and the *y*-direction implicitly and, for the second step $(n+1/2 \ \text{to}\ n+1)$, the *x*-direction is solved implicitly and the *y*-direction explicitly. The advantage of the ADI method is that the equation in each step has a simpler structure and can be solved more efficiently (e.g., with the tridiagonal matrix algorithm). Equation $(3)$ for each half-step is then given as: 

$$
\begin{equation}
\frac{T_{i,j}^{n+1/2}-T_{i,j}^n}{\Delta t/2} = 
\kappa 
    \left( 
    \frac{T_{i+1,j}^n-2T_{i,j}^n+T_{i-1,j}^n}{\Delta x^2} +
    \frac{T_{i,j+1}^{n+1/2}-2T_{i,j}^{n+1/2}+T_{i,j-1}^{n+1/2}}{\Delta y^2}
    \right)
\end{equation}
$$

$$
\begin{equation}
\frac{T_{i,j}^{n+1}-T_{i,j}^{n+1/2}}{\Delta t/2} = 
\kappa 
    \left( 
    \frac{T_{i+1,j}^{n+1}-2T_{i,j}^{n+1}+T_{i-1,j}^{n+1}}{\Delta x^2} + 
    \frac{T_{i,j+1}^{n+1/2}-2T_{i,j}^{n+1/2}+T_{i,j-1}^{n+1/2}}{\Delta y^2}
    \right)
\end{equation}
$$

&emsp; This results in two linear sets of linear system of euqations with coefficient matrices for the left and right hand side of the equations. The corresponding coefficients and right hand side of each linear system of equations needs to be adjusted according to the given boundary conditions, as shown in the CNA, for example. For more details on how this is implemented, see [*ADI.jl*](../../../src/HeatEquation/ADI.jl).

## Steady State Solution 

&emsp;So far, variable thermal parameters are only implemented in the 1-D and 2-D steady state solutions (except for the 2-D defection correction method). In steady state, one assumes that the temperature does not vary over time (i.e., $\frac{\partial T}{\partial t}=0$) and the *temperature conservation equation* simplifies to an *elliptic partial differential* equation (i.e., the *Poission equation*).  

### Poisson Solution (constant *k*)

&emsp;For constant thermal parameters the diffusive part of the *temperature conservation equation* is given by: 

$$
\begin{equation}
0 = \left( 
    \frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial z^2}
    \right) + \frac{Q}{k}.
\end{equation}
$$

&emsp;For the approximation of the spatial partial derivatives with finite difference expressions, a central finite difference is chosen and equation $(45)$ is then given as: 

$$
\begin{equation}
0 = \left( 
    \frac{T_{i+1,j} - 2T_{i,j} + T_{i-1,j}}{\Delta x^2} + \frac{T_{i,j+1} - 2T_{i,j} + T_{i,j-1}}{\Delta y^2}
    \right) + \frac{Q}{k},
    \end{equation}
$$

where *i* and *j* are the indices for the *x*- and *y*-direction, respectively. Now, one can rearrange the equation by known (*Q*, *k*) and unknown (*T*) variables, wich results in a linear system of equations in the form of: 

$$
\begin{equation} 
aT_{i-1,j} + bT_{i,j-1} - 2(a+b)T_{i,j} + bT_{i,j+1} + aT_{i+1,j}
=
-\frac{Q}{k},
\end{equation} 
$$

where $a = \frac{1}{\Delta x^2}$ and $b = \frac{1}{\Delta y^2}$. For more details on how this is implemented see [*PoissonSolvers.jl*](../../../src/HeatEquation/PoissonSolvers.jl).

&emsp;Here, one can again assume Dirichlet or Neumann boundary conditions, where we can use the temperature on the *ghost nodes* to properly implement the corresponding boundary conditions. Thereby, one needs to modify the coefficient matrix and the right hand side of the linear system of equations accordingly. For example, for a Dirichlet boundary condition on the left lateral boundary (*i = 1*) the equation changes as: 

$$
\begin{equation}
bT_{1,j-1} - (3a + 2b)T_{1,j} + bT_{1,j+1} + aT_{2,j} 
= 
-\frac{Q_{i,j}}{k_{i,j}} + 2aT_{BC,W}.
\end{equation}
$$

### Poisson Solution (variable *k*)

&emsp; For variable thermal parameters the steady-state temperature equation is given by (in 2-D): 

$$
\begin{equation} 
0 = \frac{\partial}{\partial x}\left(k_x\frac{\partial T}{\partial x}\right) + 
\frac{\partial}{\partial y}\left(k_y\frac{\partial T}{\partial y}\right) + Q_{i,j} 
\end{equation} 
$$

&emsp;To properly solve equation $(49)$, one needs to apply a conservative finite difference scheme, such that the heat flux $(q_i = k_i\frac{\partial T}{\partial i})$, and thus the thermal conductivity, is defined **between** the centroids and the temperature **on** the centroids (see Figure 1). Considering such a staggered finite difference grid, equation $(49)$ results in a linear system of equations in the form of: 

$$
\begin{equation} 
a k_{x;i,j} T_{i-1,j} + 
b k_{y;i,j} T_{i,j-1} + 
c T_{i,j} + 
b k_{y;i,j+1} T_{i,j+1} + 
a k_{x;i+1,j} T_{i+1,j} + Q_{i,j} = 0 
\end{equation}
$$

where 

$$ 
a = \frac{1}{\Delta{x^2}},
b = \frac{1}{\Delta{y^2}}, \ and \ 
c = -a\left(k_{x;i+1,j}+k_{x;i,j}\right) - b\left(k_{y;i,j+1}+k_{y;i,j}\right).
$$
For more details on how this is implemented see [*PoissonSolvers.jl*](../../../src/HeatEquation/PoissonSolvers.jl).

---------------
---------------

## Examples 

### [Backward Euler](../../../examples/HeatEquation/2D/BackwardEuler.jl)

### [Forward Euler](../../../examples/HeatEquation/2D/ForwardEuler.jl)

### [Gaussian Diffusion](../../../examples/HeatEquation/2D/Gaussian_Diffusion.jl)

### [Resolution Test Poisson Problem](../../../examples/HeatEquation/2D/Poisson_ResTest.jl)

### [Poisson Variable Parameters](../../../examples/HeatEquation/2D/Poisson_variable_k.jl) 


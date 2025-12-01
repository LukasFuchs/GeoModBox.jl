# Heat Diffusion Equation (1D)

In one dimension, the diffusive component of the *temperature evolution equation* is expressed as follows, assuming only radiogenic heat sources:

$\begin{equation}
\rho c_p\frac{\partial T}{\partial t} = -\frac{\partial q_x}{\partial x} + Q.
\end{equation}$ 

By incorporating Fourier’s law and allowing for spatially variable thermal properties, the equation becomes:

$\begin{equation}
\rho c_p\frac{\partial T}{\partial t} = \frac{\partial}{\partial x} \left( k_x \frac{\partial T}{\partial x} \right) + Q. 
\end{equation}$

Assuming constant thermal properties, Equation (2) simplifies to:

$\begin{equation}
\frac{\partial T}{\partial t} = \kappa \frac{\partial^2 T}{\partial x^2} + \frac{Q}{\rho c_p},
\end{equation}$
  
where $\kappa = k/\rho/c_p$ is the thermal diffusivity [m²/s], and $Q$ is the volumetric heat production rate [W/m³].  

Equation (3) is classified as a *parabolic partial differential equation* (PDE), which can be solved numerically given appropriate initial and boundary conditions.

# Discretization and Numerical Schemes

To solve Equation (3) numerically, the spatial domain must be discretized, assigning physical parameters to their corresponding grid locations.

> **Note**: Although thermal conductivity is currently assumed to be constant, a *conservative gridding* approach is employed to ensure physical consistency. In this scheme, temperature $T$ is defined at **cell centers** (centroids), while heat flux $q$ is defined at **cell interfaces** (vertices).

![1DDiscretization](../assets/Diff_1D_Discretization.png)

**Figure 1. 1D Discretization.** Conservative finite difference grid used to solve the 1D heat diffusion equation. Temperature is defined at centroids, while heat flux is defined at vertices. *Ghost nodes* are introduced to implement *Dirichlet* and *Neumann* boundary conditions.

The example script [Heat_1D_discretization.jl](./examples/GaussianDiffusion1D.md) demonstrates various numerical schemes for solving the heat diffusion equation, including *explicit*, *implicit*, and *Crank–Nicolson*. Below, these well-known schemes are briefly described and their respective strengths and limitations highlighted.

## Temperature Field Management

For the **Forward in Time and Centered in Space (FTCS)** explicit and implicit solvers the extended temperature field — including ghost nodes — is used to evaluate the heat diffusion equation. The old temperature field is assigned to the centroids of the extended grid to compute the new temperature.

For the additional implicit methods (Crank-Nicolson), the current temperature at the centroids is assigned to the right-hand side vector. The coefficient matrix is then assembled, and the new temperature is computed by solving the resulting linear system.

## Boundary Conditions 

Boundary conditions are implemented using ghost nodes located at $\frac{\Delta x}{2}$ outside the domain boundaries. Within the `GeoModBox.jl`, the two most common thermal boundary conditions are currently considered: Dirichlet and Neumann.

**Dirichlet Boundary Condition**

The Dirichlet condition specifies a fixed temperature at the boundary. The temperature at the left (West) and right (East) ghost nodes $T_{\textrm{G},W}$ and $T_{\textrm{G},E}$ is given by:

$\begin{equation}
T_{\textrm{G},W} = 2T_{\textrm{BC},W} - T_{1},
\end{equation}$

$\begin{equation}
T_{\textrm{G},E} = 2T_{\textrm{BC},E} - T_{\textrm{nc}},
\end{equation}$

where 
$T_{\textrm{BC},W}$ and $T_{\textrm{BC},E}$ are the prescribed boundary temperatures,
$T_1$ and $T_{\textrm{nc}}$ are the temperatures at the first and last interior centroids, and
$\textrm{nc}$ is the number of internal centroids.

**Neumann Boundary Condition**

The Neumann condition specifies a fixed gradient (e.g., a heat flux or temperature gradient) across the boundary. The ghost node temperatures are defined as:

$\begin{equation}
T_{\textrm{G},W} = T_{1} - c_{W} \Delta{x},
\end{equation}$

$\begin{equation}
T_{\textrm{G},E} = T_{\textrm{nc}} + c_{E} \Delta{x},
\end{equation}$

with:  

$\begin{equation}
\left. c_{W} = \frac{\partial{T}}{\partial{x}} \right\vert_{W},\ \textrm{and}\ \left. c_{E} = \frac{\partial{T}}{\partial{x}} \right\vert_{E}, 
\end{equation}$

where $c_W$ and $c_E$ are the prescribed temperature gradients across the west and east boundaries, respectively. These ghost node values can be substituted into the discretized heat diffusion equation for the internal centroids adjacent to the boundaries. Using ghost nodes consistently enforces the boundary conditions at each time step and preserves the second-order accuracy treatment of the three-point stencil (e.g. Duretz et al., 2011).

## Explicit Finite Difference Scheme (FTCS; Forward Euler)

A fundamental and intuitive approach solving the 1D heat diffusion equation is the **FTCS** scheme, implemented in an **explicit** manner.

This method approximates the continuous PDE on a discrete grid and converges to the analytical solution as the spatial ($\Delta x$) and temporal ($\Delta t$) resolutions are refined. Its main advantages are **simplicity** and **computational efficiency**.

However, the FTCS scheme is **conditionally stable**. Its stability is governed by the *heat diffusion stability criterion*, which can be derived via *Von Neumann* analysis. This assesses how numerical perturbations grow or decay over time. 

In the following, we implement the FTCS scheme by combining a three-point stencil in space with a two-point stencil in time. 

For a uniform, one-dimensional grid, the stability condition is:

$\begin{equation}
\Delta t < \frac{\Delta{x^2}}{2 \kappa}.
\end{equation}$ 

Consequently, the maximum allowable time step is constrained by the spatial resolution.

Discretizing Equation (3) with the FTCS scheme gives:

$\begin{equation}
\frac{T_{i}^{n+1} - T_{i}^{n} }{\Delta t} = \kappa \frac{T_{i-1}^{n} - 2T_{i}^{n} + T_{i+1}^{n}}{\Delta{x^2}} + \frac{Q_{i}^n}{\rho c_p},
\end{equation}$ 

where 
$i$ is the spatial grid index,
$n$ is the time step index,
$\Delta x$ is the grid spacing, and
$\Delta t$ is the time step.

Solving for $T_i^{n+1}$ gives 

$\begin{equation}
T_{i}^{n+1} = T_{i}^{n} + a \left(T_{i-1}^{n} - 2T_{i}^{n} + T_{i+1}^{n} \right) + \frac{Q_{i}^n \Delta t}{\rho c_p}, 
\end{equation}$

where 

$\begin{equation} 
a = \frac{\kappa \Delta t}{\Delta x^2}.
\end{equation}$

Equation (11) is solved for all interior nodes at each time step, assuming initial and boundary conditions are specified. The boundary conditions are implemented using the temperature values calculated for the ghost nodes. For implementation details, see the [source code](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/HeatEquation/1Dsolvers.jl).

## Implicit Scheme (Backward Euler)

The fully implicit finite difference scheme, also known as the **Backward Euler** method, is **unconditionally stable**, allowing time steps larger than those permitted by the CFL criterion.

> **Note**: While the implicit method is unconditionally stable, very large time steps may still yield to larger inaccuracies with respect to an analytical solution, particularly for resolving small-scale thermal gradients (THIS STATEMENT IS NOT CLEAR: which scale? spatial scale? does running with poor time resolution strongly afftec spatial resolution?).

In 1D and based on a three-point stencil, the discretized heat diffusion equation becomes:

$\begin{equation}
\frac{T_{i}^{n+1}-T_{i}^n}{\Delta t} = \kappa \frac{T_{i-1}^{n+1}-2T_{i}^{n+1}+T_{i+1}^{n+1}}{\Delta{x^2}} + \frac{Q^n}{\rho c_p},
\end{equation}$

where $i$ is the spatial index, $n$ is the time step index, $\Delta t$ is the time step, and $\Delta x$ is the spatial grid spacing. Rearranging the equation into known (right-hand side) and unknown (left-hand side) terms yields a system of equations considering all internal centroids:

$\begin{equation}
-a T_{i-1}^{n+1} + \left(2a + b \right) T_{i}^{n+1} - a T_{i+1}^{n+1} = b T_{i}^n + \frac{Q^n}{\rho c_p},
\end{equation}$

with 

$\begin{equation}\begin{split}
a & = \dfrac{\kappa}{\Delta x^2} \text{, and} \\
b & = \dfrac{1}{\Delta t}.
\end{split}\end{equation}$

These equations are a tridiagonal system of equations in the form of:

$\begin{equation}
\mathbf{K} \cdot \bm{x} = \bm{b}
\end{equation}$

where $\mathbf{K}$ is the coefficient matrix (with three non-zero diagonals), $\bm{x}$ is the unknown solution vector, that is the temperature at time step $n+1$, and $\bm{b}$ is the known right-hand side.

### General Solution

A general approach solving this system of equations is the **defection correction**. The heat diffusion equation is reformulated by introducing a residual term $\bm{r}$, which quantifies the deviation from the true solution and can be reduced iteratively to improve accuracy through successive correction steps. In implicit form, Equation (16) can be rewritten as:

$\begin{equation}
\mathbf{K} \cdot \bm{x} - \bm{b} = \bm{r}, 
\end{equation}$

where $\bm{r}$ is the residual (or defect). The coefficients of the matrix are the same as derived above, but can now also be calculated by the Jacobian of $\bm{r}$ via: 

$\begin{equation}
{K_{ij}}=\frac{\partial{{r}_i}}{\partial{{x}_j}}.
\end{equation}$

We can solve for the unknown vector $\bm{x}$ by assuming an initial temperature guess, calculating its residual, and performing a correction of the temperature. For non-linear problems, the process is repeated iteratively until the residual is sufficiently small. In the linear case, one iteration yields the exact solution.

Given an initial temperature guess $\bm{T}^k$, the initial residual is:

$\begin{equation}
\bm{r}^k = \mathbf{K} \cdot \bm{T}^k - \bm{b}.
\end{equation}$

To reduce the residual, a correction term $\delta \bm{T}$ is defined such that:

$\begin{equation}
0 = \bm{K} \left(\bm{T}^k + \delta \bm{T} \right) - \bm{b} = \bm{K} \bm{T}^k - \bm{b} + \bm{K} \delta \bm{T} = \bm{r}^k + \bm{K} \delta \bm{T}.
\end{equation}$

Rearranging gives: 

$\begin{equation}
\bm{r}^k = -\mathbf{K} \delta \bm{T}, 
\end{equation}$

and hence: 

$\begin{equation}
\delta \bm{T} = -\mathbf{K}^{-1} \bm{r}^k. 
\end{equation}$

Thus, the updated solution becomes: 

$\begin{equation}
\bm{T}^{k+1} = \bm{T}^k + \delta \bm{T},
\end{equation}$

where $\bm{T}^{k+1}$ is the updated temperature after one iteration step. 

Within `GeoModBox.jl` the residual $\bm{r}$ is calculated on the internal centroids using the extended temperature field including the ghost nodes of the current time step as an initial guess: 

$\begin{equation}
\frac{\partial{T_{\textrm{ext}}}}{\partial{t}} - \kappa \frac{\partial^2{T_{\textrm{ext}}}}{\partial{x^2}} - \frac{Q}{\rho c_p} = r,  
\end{equation}$

and in discretized, implicit finite-difference form: 

$\begin{equation}
\frac{T_{\textrm{ext},i}^{n+1} - T_{\textrm{ext},i}^{n}}{\Delta{t}} - \kappa \frac{T_{\textrm{ext},i-1}^{n+1} - 2 T_{\textrm{ext},i}^{n+1} + T_{\textrm{ext},i+1}^{n+1}}{\Delta{x^2}} - \frac{Q^n}{\rho c_p}= r_j, 
\end{equation}$

where $i=2:(nc+1)$ is the index of the internal centroid of the extended temperature field and $j=1:nc$ is the index of the internal centroids of the regular field. Rewriting Equation (25) and substituting the coefficients using Equation (15) results in: 

$\begin{equation}
-aT_{\textrm{ext},i-1}^{n+1}
+\left(2a+b\right)T_{\textrm{ext},i}^{n+1}
-aT_{\textrm{ext},i+1}^{n+1}
-bT_{\textrm{ext},i}^n
-\frac{Q^n}{\rho c_p}=r_j,
\end{equation}$

which corresponds to the matrix form of Equation (17), where $T_{\textrm{ext},i}^{n+1}$ is the unknown vector $x$,$-bT_{\textrm{ext},i}^n -\frac{Q^n}{\rho c_p}$ is the known vector $b$, and $-a$ and $2a+b$ are the coefficients of the non-zero diagonals of the coefficient matrix. With the residual vector $r$ and the coefficient matrix $\bm{K}$ one can calculate the correction term for the temperature via Equation (22). 

### Boundary Condition 

The boundary conditions are implemented again using the temperatures values at the ghost nodes (see Equations (4)-(7)). To maintain symmetry in the coefficient matrix, however, the matrix coefficients must be modified for the nodes adjacent to the boundaries.

The equations for the internal centroids adjacent to the boundary are then given by: 

**Dirichlet Boundary Condition**

**West Boundary**

$\begin{equation}
-aT_{\textrm{G},W}
+\left(2a+b\right)T_{\textrm{ext},2}^{n+1}
-aT_{\textrm{ext},3}^{n+1}
-bT_{\textrm{ext},2}^n
-\frac{Q^n}{\rho c_p}=r_1,
\end{equation}$

$\begin{equation}
-a\left(2T_{\textrm{BC},W}-T_{\textrm{ext},2}^{n+1}\right)
+\left(2a+b\right)T_{\textrm{ext},2}^{n+1}
-aT_{\textrm{ext},3}^{n+1}
-bT_{\textrm{ext},2}^n
-\frac{Q^n}{\rho c_p}=r_1,
\end{equation}$

$\begin{equation}
\left(3a+b\right)T_{\textrm{ext},2}^{n+1}
-aT_{\textrm{ext},3}^{n+1}
-bT_{\textrm{ext},2}^n
-2aT_{\textrm{BC},W}
-\frac{Q^n}{\rho c_p}=r_1,
\end{equation}$

**East Boundary**

$\begin{equation}
-aT_{\textrm{ext},nc}^{n+1}
+\left(2a+b\right)T_{\textrm{ext},nc+1}^{n+1}
-aT_{\textrm{G},E}
-bT_{\textrm{ext},nc+1}^n
-\frac{Q^n}{\rho c_p}=r_{\textrm{nc}},
\end{equation}$

$\begin{equation}
-aT_{\textrm{ext},nc}^{n+1}
+\left(2a+b\right)T_{\textrm{ext},nc+1}^{n+1}
-a\left(2T_{\textrm{BC},E}-T_{\textrm{ext},nc+1}^{n+1}\right)
-bT_{\textrm{ext},nc+1}^n
-\frac{Q^n}{\rho c_p}=r_{\textrm{nc}},
\end{equation}$

$\begin{equation}
-aT_{\textrm{ext},nc}^{n+1}
+\left(3a+b\right)T_{\textrm{ext},nc+1}^{n+1}
-bT_{\textrm{ext},nc+1}^n
-2aT_{\textrm{BC},E}
-\frac{Q^n}{\rho c_p}=r_{\textrm{nc}},
\end{equation}$

**Neumann Boundary Conditions**

**West Boundary**

$\begin{equation}
-aT_{\textrm{G},W}
+\left(2a+b\right)T_{\textrm{ext},2}^{n+1}
-aT_{\textrm{ext},3}^{n+1}
-bT_{\textrm{ext},2}^n
-\frac{Q^n}{\rho c_p}=r_1, 
\end{equation}$

$\begin{equation}
-a\left(T_{\textrm{ext},2}^{n+1} - c_{W} \Delta{x}\right)
+\left(2a+b\right)T_{\textrm{ext},2}^{n+1}
-aT_{\textrm{ext},3}^{n+1}
-bT_{\textrm{ext},2}^n
-\frac{Q^n}{\rho c_p}=r_1, 
\end{equation}$

$\begin{equation}
+\left(a+b\right)T_{\textrm{ext},2}^{n+1}
-aT_{\textrm{ext},3}^{n+1}
-bT_{\textrm{ext},2}^n
+a c_{W} \Delta{x}
-\frac{Q^n}{\rho c_p}=r_1, 
\end{equation}$

**East Boundary**

$\begin{equation}
-aT_{\textrm{ext},nc}^{n+1}
+\left(2a+b\right)T_{\textrm{ext},nc+1}^{n+1}
-aT_{\textrm{G},E}
-bT_{\textrm{ext},nc+1}^n
-\frac{Q^n}{\rho c_p}=r_{\textrm{nc}},
\end{equation}$

$\begin{equation}
-aT_{\textrm{ext},nc}^{n+1}
+\left(2a+b\right)T_{\textrm{ext},nc+1}^{n+1}
-a\left(T_{\textrm{ext},nc+1}^{n+1} + c_{E} \Delta{x}\right)
-bT_{\textrm{ext},nc+1}^n
-\frac{Q^n}{\rho c_p}=r_{\textrm{nc}},
\end{equation}$

$\begin{equation}
-aT_{\textrm{ext},nc}^{n+1}
+\left(a+b\right)T_{\textrm{ext},nc+1}^{n+1}
-bT_{\textrm{ext},nc+1}^n
-ac_{E} \Delta{x}
-\frac{Q^n}{\rho c_p}=r_{\textrm{nc}},
\end{equation}$

where 
$T_{\textrm{BC},W}$, $T_{\textrm{BC},E}$ are the prescribed boundary temperatures,
$c_W$, $c_E$ are the prescribed temperature gradients at the west and east boundaries, and
$nc$ is the number of internal centroids. 

These adjustments ensure that the boundary conditions are enforced consistently while preserving the symmetry of the implicit solver.

### Special Case - A Linear Problem 

If the problem is linear and the exact solution is reached within one single iteration step, the system of equations is reduced to Equation (16). Thus, one can solve the system of equations directly via a *left matrix division*: 

$\begin{equation}
\bf{x} = \mathbf{K}^{-1} \bm{b}. 
\end{equation}$

The coefficient matrix remains the same, even for the given boundary conditions. However, the right-hand side needs to be updated accordingly (simply setting $\bm{r}=0$ and adding the known parameters to the right-hand side of the equations). Thus, the equations for the internal centroids adjacent to the boundaries are defined as: 

**Dirichlet Boundary Condition**

**West boundary**

$\begin{equation}
\left(3 a + b\right) T_{1}^{n+1} - a T_{2}^{n+1} = b T_{1}^{n} + 2 a T_{\textrm{BC},W},
\end{equation}$

**East boundary**

$\begin{equation}
-a T_{\textrm{nc}-1}^{n+1} + \left(3 a + b\right) T_{\textrm{nc}}^{n+1}  = b T_{\textrm{nc}}^{n} + 2 a T_{\textrm{BC},E}, 
\end{equation}$

**Neumann Boundary Condition**

**West boundary**

$\begin{equation}
\left(a + b\right) T_{1}^{n+1} - a T_{2}^{n+1} = b T_{1}^{n} - a c_{W} \Delta{x},
\end{equation}$

**East boundary**

$\begin{equation}
-a T_{\textrm{nc}-1}^{n+1} + \left(a + b\right) T_{\textrm{nc}}^{n+1}  = b T_{\textrm{nc}}^{n} + a c_{E} \Delta{x}. 
\end{equation}$

For the sake of simplicity and since most of the problems within `GeoModBox.jl` so far are of linear nature, the *special case* formulation is applied to the additional methods. 

## Crank-Nicolson approach (CNA)

The fully implicit FTCS method is unconditionally stable but only first-order accurate in time. To improve temporal accuracy while retaining stability, the **Crank-Nicolson scheme** can be used. This method employs a time-centered (implicit) discretization and is second-order accurate in time.

In one dimension, the Crank-Nicolson discretization of the heat diffusion equation becomes:

$\begin{equation}
\frac{T_{i}^{n+1} - T_{i}^{n}}{\Delta t} = \frac{\kappa}{2}\frac{(T_{i-1}^{n+1}-2T_{i}^{n+1}+T_{i+1}^{n+1})+(T_{i-1}^{n}-2T_{i}^{n}+T_{i+1}^{n})}{\Delta{x^2}}. 
\end{equation}$

Rearranging into known and unknown terms yields a linear system of the form:

$\begin{equation}
-aT_{i-1}^{n+1} + \left(b+2a\right)T_{i}^{n+1} - a T_{i+1}^{n+1} = aT_{i-1}^{n} + \left(b-2a\right)T_{i}^{n} + a T_{i+1}^{n},
\end{equation}$

where:

$\begin{equation}\begin{split}
a & = \frac{\kappa}{2\Delta{x^2}} \text{, and} \\
b & = \frac{1}{\Delta{t}}.
\end{split}\end{equation}$

### Boundary Conditions

To obtain a symmetric coefficient matrix, both the matrix and the right-hand side vector must be modified at the boundaries. The equations for centroids adjacent to the boundaries are:

**Dirichlet Boundary Conditions**

**West boundary**

$\begin{equation}
\left(b + 3 a \right) T_{1}^{n+1} - a T_{2}^{n+1} = \left( b - 3 a \right) T_{1}^{n} + a T_{2}^{n} + 4 a T_{\textrm{BC},W}
\end{equation}$

**East boundary**

$\begin{equation}
-a T_{\textrm{nc}-1}^{n+1} + \left(b + 3 a \right) T_{\textrm{nc}}^{n+1} = a T_{nc-1}^{n} + \left( b - 3 a \right) T_{\textrm{nc}}^{n} + 4 a T_{\textrm{BC},E}
\end{equation}$

**Neumann Boundary Conditions**

**West boundary**

$\begin{equation}
\left(b+a\right)T_{1}^{n+1} - a T_{2}^{n+1} = \left(b-a\right)T_{1}^{n} + a T_{2} - 2ac_{W} \Delta{x}
\end{equation}$

**East boundary**

$\begin{equation}
-a T_{\textrm{nc}-1}^{n+1} + \left(b+a\right)T_{\textrm{nc}}^{n+1}  = a T_{\textrm{nc}-1}^{n} + \left(b-a\right)T_{\textrm{nc}}^{n} + 2ac_{E} \Delta{x}
\end{equation}$

The resulting coefficient matrix remains tridiagonal, preserving computational efficiency. However, memory requirements increase with finer spatial resolution due to the larger size of the linear system, making this method more memory intensive.

For implementation details, refer to the [source code](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/HeatEquation/1Dsolvers.jl).

## Summary

While the **explicit FTCS scheme** is simple and efficient for small time steps, **implicit methods** like Backward Euler and Crank-Nicolson are preferred for their unconditional stability. The **defection correction** provides a flexible framework for both linear and nonlinear problems, allowing for iterative refinement when needed. The Crank-Nicolson scheme further improves accuracy with its second-order time discretization. 

---

# Variable Thermal Parameters 

To solve the 1D head diffusion equation with **spatially variable thermal properties**, a conservative finite difference scheme is employed. In this formulation, temperature is defined at centroids, while heat flux and thermal conductivity are defined at vertices (see Figure 1).

>Note: Currently, `GeoModBox.jl` only provides an explicit solver for the 1-D solution of the heat diffusion equation including variable thermal parameters. Thus, only the explicit scheme is described in the following. 

The governing equation is:

$\begin{equation}
\rho c_{p} \frac{\partial{T}}{\partial{t}} = \frac{\partial{}}{\partial{y}}\left(k \frac{\partial{T}}{\partial{y}}\right) + \rho H,
\end{equation}$ 

where
$\rho$ is the density [kg/m³],
$c_p$ is the specific heat capacity [J/kg/K],
$T$ is the temperature [K],
$t$ is the time [s],
$k$ is the thermal conductivity [W/m/K], 
$H$ is the internal heat generation rate per unit mass [W/kg], and
$y$ is the vertical coordinate (depth) [m]

### Discretization

In a conservative scheme, the vertical conductive heat flux $q_y$ is defined on vertices, as:

$\begin{equation}
\left. q_{\textrm{y},m} = -k_m \frac{\partial T}{\partial y}\right\vert_{m},\ \textrm{for}\ m = 1:nv, 
\end{equation}$

where $nv$ is the number of *vertices*.

### Explicit Finite Difference Scheme

Using the above discretization, the heat diffusion equation at each centroid is computed from:

$\begin{equation}
\rho_j c_{p,j} \frac{T_{j}^{n+1} - T_{j}^{n}}{\Delta{t}} = -\frac{q_{y,j+1}^{n} - q_{y,j}^{n} }{\Delta{y}} + \rho_j H_j,\ \textrm{for}\ j = 1:nc, 
\end{equation}$

where 
$T_j$ is evaluated at centroids,
$q_y$ and $k$ are evaluated at vertices,
$\Delta t$ is the time step, and
$\Delta y$ is the spatial grid resolution.

Substituting the expression for $q_y$ gives:

$\begin{equation}
\rho_j c_{p,j} \frac{T_{j}^{n+1} - T_{j}^{n}}{\Delta{t}} = \frac{ k_{j+1} \frac{T_{j+1}^{n} - T_{j}^{n}}{\Delta{y}} - k_{j} \frac{T_{j}^{n} - T_{j-1}^{n}}{\Delta{y}} }{\Delta{y}} + \rho_j H_j.
\end{equation}$

Rewriting this into an explicit update formula for $T_j^{n+1}$:

$\begin{equation}
T_{j}^{n+1} = ak_{j}T_{j-1}^{n} + \left(1-a\left(k_{j+1}+k_{j}\right)\right)T_{j}^{n} + ak_{j+1}T_{j+1}^{n} + \frac{H_j\Delta{t}}{c_{p,j}},
\end{equation}$

where:

$\begin{equation}
a = \frac{\Delta{t}}{\Delta{y^2} \rho c_{p_{j}}}.
\end{equation}$

### Boundary Conditions

For centroids adjacent to the boundaries, ghost nodes are used to evaluate the temperature gradient consistently with the chosen thermal boundary condition (Dirichlet or Neumann). These ghost node values are computed according to equations (7)–(10).

---

For implementation details, refer to the [source code](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/HeatEquation/1Dsolvers.jl).

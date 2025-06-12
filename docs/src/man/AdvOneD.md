# Advection Equation (1D)

In one dimension, the advection for the temperature conservation equation, for example, is given as follows

$\begin{equation}
\frac{\partial{T}}{\partial{t}} = -v_x \left(\frac{\partial{T}}{\partial{x}}\right),
\end{equation}$

where $T$ is the temperature [K], $t$ is the time [s], and $v_x$ is the velocity in $x$-direction. 

... *Courant time step* ...

# Discretization Schemes

## Forward in Time and Centered in Space (FTCS)

Let's first use the seemingly easiest approach to solve the advection equation. Approximating the partial derivative operators using a *FTCS* scheme results in 

$\begin{equation}
T_i^{n+1} = T_i^n - v_x \Delta{t}\frac{T_{i+1}^n - T_{i-1}^n}{2\Delta{x}},
\end{equation}$

where $\Delta{t}$ and $\Delta{x}$ are the time step and grid resolution, respectively, and $i$ and $n$ are the index in the $x$-direction and the current time step. 

Unfortunately, this scheme is always unstable for the advection equation, as can be shown by a *Von Neumann or Hirt's stability analysis* [^1]. Besides these, one can show, that the central difference a $i$ always results in an amplification of the information (here the temperature) at $i$ at the next time step $n+1$. Thus, it is continually increasing and *unconditionally unstabel*. 

A way out of this dilema, is to only consider the information from upstream. 

## Upwind

The upwind scheme only uses one-sided, spatial finite difference approximations always taken in the upstream (upwind) direction. That is, one needs to take the actual velocity into account. The discretized form of the advection equation is then given as 

$\begin{equation}
\frac{T_{i}^{n+1}-T_{i}^n}{\Delta{t}} = -v_{x,i}
\begin{cases}
\frac{T_{i}^{n}-T_{i-1}^{n}}{\Delta{x}} &\text{if } v_{x,i} \gt 0\\
\frac{T_{i+1}^{n}-T_{i}^{n}}{\Delta{x}}&\text{if } v_{x,i} \lt 0 
\end{cases}.
\end{equation}$

The scheme is stable as long as the *Courant* criterion is satisfied, but it has a strong numerical diffusion which depends on the grid size. Using a Taylor series one can show, that, in 1-D and for a constant velocity, no numerical diffusion occurs if the time step is exactly the *Courant* criterion. The method is unstable, if the *Courant* criterion is violated. 

## LAX

Modification of the *FTCS* scheme by damping the instability

- Lax Method -> stable but strong numerical diffusion for $\alpha \ne 1$

$\begin{equation}
\frac{T_{i}^{n+1}-\left(T_{i+1}^{n}+T_{i-1}^{n}\right)/2}{\Delta{t}}=-v_x\frac{T_{i+1}^{n}-T_{i-1}^{n}}{2\Delta{x}}
\end{equation}$

## Staggered Leapfrog

- Staggered Leap Frog -> Higher Fehlerordnung in time and space, no numerical diffusion 

$\begin{equation}
\frac{T_{i}^{n+1}-T_{i}^{n-1}}{2\Delta{t}}=-v_x\frac{T_{i+1}^{n}-T_{i-1}^{n}}{2\Delta{x}}
\end{equation}$

## Semi-Lagrangian Method

- Semi-lagrangian

$\begin{equation}
\end{equation}$

## Passive tracers
- Tracer

$\begin{equation}
\end{equation}$

## Exercise

All the above listed methods are part of the 1-D advection exercise. 

- [1-D Gaussian or block anomaly advection](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/06_1D_Advection.ipynb)  

# References

[^1]: Spiegelman, M. (2004). Myths and methods in modeling. Columbia University Course Lecture Notes, available online at http://www. ldeo. columbia. edu/~ mspieg/mmm/course. pdf, accessed, 6, 2006.
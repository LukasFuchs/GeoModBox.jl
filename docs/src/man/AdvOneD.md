# Advection Equation (1D)

In one dimension, the advection for the temperature conservation equation, for example, is given as follows

$\begin{equation}
\frac{\partial{T}}{\partial{t}} = -v_x \left(\frac{\partial{T}}{\partial{x}}\right),
\end{equation}$

where $T$ is the temperature [K], $t$ is the time [s], and $v_x$ is the velocity in $x$-direction. 

# Discretization Schemes

## Forward in Time and Centered in Space (FTCS)

Let's first use the seemingly easiest approach to solve the advection equation. Approximating the partial derivative operators using a *FTCS* scheme results in 

$\begin{equation}
T_i^{n+1} = T_i^n - v_x \Delta{t}\frac{T_{i+1}^n - T_{i-1}^n}{2\Delta{x}},
\end{equation}$

where $\Delta{t}$ and $\Delta{x}$ are the time step and grid resolution, respectively, and $i$ and $n$ are the index in the $x$-direction and the current time step. 

Unfortunately, this scheme is always unstable for the advection equation! 

... *Why? tba* ...

## Upwind

- Upwind -> Courtant criterium, numerische diffusion, exact solution for \alpha = 1
    -> error estimaiton of each FD term, taylor expansion 
    -> \alpha = 1 no numerical diffusion, resolution dependent

$\begin{equation}
\frac{T_{i}^{n+1}-T_{i}^n}{\Delta{t}} = -v_{x,i}
\begin{cases}
\frac{T_{i}^{n}-T_{i-1}^{n}}{\Delta{x}} &\text{if } v_{x,i} \gt 0\\
\frac{T_{i+1}^{n}-T_{i}^{n}}{\Delta{x}}&\text{if } v_{x,i} \lt 0 
\end{cases}.
\end{equation}$

## LAX

- Lax Method -> stable but strong numerical diffusion

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

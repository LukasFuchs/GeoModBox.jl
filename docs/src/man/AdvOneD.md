# Advection Equation (1D)

In one dimension, the advection equation for the temperature conservation, for example, is given as follows

$\begin{equation}
\frac{\partial{T}}{\partial{t}} = -v_x \left(\frac{\partial{T}}{\partial{x}}\right),
\end{equation}$

where $T$ is the temperature [K], $t$ is the time [s], and $v_x$ is the velocity in the $x$-direction. 

# Discretization Schemes

## Forward in Time and Centered in Space (FTCS)

Let's begin with the seemingly simplest approach. Approximating the partial derivatives using a *FTCS* scheme results in

$\begin{equation}
\frac{T_{i}^{n+1}-T_{i}^{n}}{\Delta{t}} = -v_x\left(\frac{T_{i+1}^{n}-T_{i-1}^{n}}{2\Delta{x}}\right),
\end{equation}$

where $\Delta{t}$ and $\Delta{x}$ are the time step and grid resolution, respectively, and $i$ and $n$ denote the spatial and temporal indices. This scheme is first-order accurate in time and second-order accurate in space. 

Rearranging gives the solution for the temperature at the next time step:

$\begin{equation}
T_i^{n+1} = T_i^n - v_x \Delta{t}\frac{T_{i+1}^n - T_{i-1}^n}{2\Delta{x}}.
\end{equation}$

The right-hand side can be simplified using the so-called *Courant number*:

$\begin{equation}
\alpha = \frac{v_x\Delta{t}}{\Delta{x}},
\end{equation}$

which represents the number of grid points traversed in a single time step. 

Unfortunately, this scheme is unconditionally unstable for the advection equation, as shown by a *Von Neumann* or *Hirt's stability analysis*. Furthermore, the central difference at $i$ causes amplification of the variable (here, temperature) at each subsequent time step. Hence, the solution continually grows and is unstable.

## Lax-Friedrichs method

One way to suppress the instability of the FTCS scheme is the *Lax-Friedrichs* method. This replaces the term $T_{i}^{n}$ with its spatial average at the same time level, resulting in:

$\begin{equation}
\frac{T_{i}^{n+1}-\left(T_{i+1}^{n}+T_{i-1}^{n}\right)/2}{\Delta{t}}=-v_x\frac{T_{i+1}^{n}-T_{i-1}^{n}}{2\Delta{x}}.
\end{equation}$

Rearanging gives:

$\begin{equation}
T_{i}^{n+1} = \frac{1}{2}\left(T_{i+1}^{n}+T_{i-1}^{n}\right)-
\frac{v_x \Delta{t}}{2\Delta{x}} \left(T_{i+1}^{n}-T_{i-1}^{n}\right).
\end{equation}$

This method is stable for $\alpha < 1$ but introduces significant numerical diffusion.

## Upwind

Another approach is to consider only upstream information. The *upwind* scheme uses one-sided finite differences, always taken in the upstream direction. This results in a scheme that is first-order accurate in both space and time. The discretized advection equation becomes:

$\begin{equation}
\frac{T_{i}^{n+1}-T_{i}^n}{\Delta{t}} = -v_{x,i}
\begin{cases}
\frac{T_{i}^{n}-T_{i-1}^{n}}{\Delta{x}} &\text{if } v_{x,i} \gt 0\\
\frac{T_{i+1}^{n}-T_{i}^{n}}{\Delta{x}}&\text{if } v_{x,i} \lt 0 
\end{cases}.
\end{equation}$

The scheme is stable if the *Courant criterion* is satisfied ($\alpha \le 1$), but numerical diffusion remains, which depends on the grid size. A Taylor series expansion shows that, in 1D with constant velocity, the scheme becomes non-diffusive if the time step exactly satisfies the Courant criterion. The method becomes unstable if this criterion is violated.

## Staggered Leapfrog

All previously discussed explicit schemes are only first-order accurate in time and second-order in space (except upwind, which is first-order in both). To match the temporal and spatial accuracy without choosing a very small time step, one may use the *staggered leapfrog* scheme:

$\begin{equation}
\frac{T_{i}^{n+1}-T_{i}^{n-1}}{2\Delta{t}}=-v_x\frac{T_{i+1}^{n}-T_{i-1}^{n}}{2\Delta{x}}.
\end{equation}$

This method avoids numerical diffusion, but becomes increasingly unstable when strong gradients in the advected field are present.

## Semi-Lagrangian Method

The methods discussed above each have drawbacks. The *semi-Lagrangian* method addresses several of them: it is stable, does not suffer from numerical diffusion, and is not constrained by the Courant criterion. It is related to tracer-based advection schemes and solves ODEs rather than using traditional finite differences. While not inherently conservative and subject to minor interpolation errors, it offers promising accuracy and efficiency.

The central idea is to trace an advected particle backward in time to its origin and interpolate the corresponding value from the Eulerian grid.

In 1D, assuming constant velocity in time and space, the procedure is:

**1. Calculate the initial position** 

The initial position $X_i$ of a particle landing on the Eulerian grid point $x_i$ at time $t_{n+1}$ is:

$\begin{equation}
X_i=x_i-\Delta{t}\cdot v_x\left(t_{n+1},x_i\right),
\end{equation}$

where $x_i$ is the coordinate of the Eulerian grid point $i$, $\Delta{t}$ is the time step, $v_x$ the velocity in $x$-direction, and $t_{n+1}$ is the time at the new time step. 

**2. Interpolate the temperature**

Interpolate the temperature at $t_n$ from the surrounding Eulerian grid points onto the position $X_i$, e.g., using `cubic_spline_interpolation()`.

**3. Update the Temperature field**

Assuming the temperature at the grid point at $t_{n+1}$ equals the interpolated value at $X_i$ at time $t_n$:

$\begin{equation}
T\left(t_{n+1},x_i\right) = T\left(t_n,X_i\right)
\end{equation}$

## Passive tracers

Using passive tracers represents a fully Lagrangian method. Initially distributed tracers (or markers) are transported by the prescribed velocity field. Tracers can carry various attributes, and if these influence the model’s rheology or dynamics, they are considered *active tracers*, requiring velocity correction.

To transport a property with tracers:

**1. Define tracers** $\vec{x}_p$ with initial positions $\vec{x}_p\left(t=0\right)$ and initial property values $f\left(\vec{x}_p\left(t=0\right)\right)$.

**2. Compute flow paths** by solving the ODE of particle motion, for instance using Forward Euler or Runge-Kutta integration.

In 1D, the path equation is:

$\begin{equation}
\frac{dx_p}{dt}=v_x \left(x_p,t\right),
\end{equation}$

where $x_p$ is the $x$-coordinate of the tracer. 

**Forward Euler**

The flow path ODE is approximated as:

$\begin{equation}
\frac{x_p^{n+1}-x_p^n}{\Delta{t}} = v_x(x_p^n). 
\end{equation}$

 Solving for the next position:

$\begin{equation}
x_p^{n+1} = x_p^n + \Delta{t}\cdot v_x(x_p^n). 
\end{equation}$

While simple, this method suffers from inaccuracy for large $\Delta{t}$ and $v_x$.
    
**Runge-Kutta 4-th order**

A more accurate method is the *4th-order Runge-Kutta*. In 1D, the next position is:

$\begin{equation}
x_p^{n+1} = x_p^n + \frac{1}{6}k_1 + \frac{1}{3}k_2 + \frac{1}{3}k_3 + \frac{1}{6}k_4,
\end{equation}$

where:

$\begin{equation}
\begin{split}
k_1 & = \Delta{t} \cdot v_x(t^n,x_p^n) \\
k_2 & = \Delta{t} \cdot v_x(t^n+\Delta{t}/2,x_p^n+k_1/2) \\
k_3 & = \Delta{t} \cdot v_x(t^n+\Delta{t}/2,x_p^n + k_2/2) \\
k_4 & = \Delta{t} \cdot v_x(t^n+\Delta{t},x_p^n+k_3) \\
\end{split}
\end{equation}$

**3. Interpolate grid values** of $f(x,t)$ from the tracer positions $\vec{x}_p$, e.g., using bilinear interpolation.

Despite the advantages, care is required. Interpolation between grid and tracer data can cause smoothing and numerical diffusion, particularly in regions with sharp gradients. Additionally, clustering or depletion of tracers can introduce further errors and may require adaptive insertion of new tracers in under-sampled regions.

## Exercise

All the above methods are implemented in the 1D advection exercise:

- [1-D Gaussian or block anomaly advection](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/06_1D_Advection.ipynb)  

# References

Spiegelman, M. (2004). Myths and methods in modeling. Columbia University Course Lecture Notes, available online at http://www. ldeo. columbia. edu/~ mspieg/mmm/course. pdf, accessed, 6, 2006.

W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, Numerical Recipes 1986, (Cambridge Univ. Press, Cambridge, 1986).
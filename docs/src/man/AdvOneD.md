# Advection Equation (1D)

In one dimension, the advection equation for the temperature conservation, for example, is given as follows

$\begin{equation}
\frac{\partial{T}}{\partial{t}} = -v_x \left(\frac{\partial{T}}{\partial{x}}\right),
\end{equation}$

where $T$ is the temperature [K], $t$ is the time [s], and $v_x$ is the velocity in the $x$-direction. 

# Discretization Schemes

The global indexing of the central reference point $I^C$ used in the advection equation follows the indexing as described in the [general solution section](./GESolution.md). The indices of the adjacent points are then defined by:

$\begin{equation}\begin{split}
I^\textrm{W} & = I^\textrm{C} - 1,\\\   
I^\textrm{E} & = I^\textrm{C} + 1,
\end{split}\end{equation}$

where $I$ is the equation number, which corresponds to the local index $i$ and the central position of the three-point stencil ($C$), and $I^\textrm{W}$, and $I^\textrm{E}$, are the points West and East of it. These are the indices of the three-point stencil used in the discretized FD equations below.

## Forward in Time and Centered in Space (FTCS)

Let's begin with the seemingly simplest approach. Approximating the partial derivatives using a *FTCS* scheme results in

$\begin{equation}
\frac{T_{I^C}^{n+1}-T_{I^C}^{n}}{\Delta{t}} = -v_x\left(\frac{T_{I^E}^{n}-T_{I^W}^{n}}{2\Delta{x}}\right),
\end{equation}$

where $\Delta{t}$ and $\Delta{x}$ are the time step and grid resolution, respectively, $I^C$ is the central reference point, and $n$ denotes the temporal index. This scheme is first-order accurate in time and second-order accurate in space. 

Rearranging gives the solution for the temperature at the next time step:

$\begin{equation}
T_{I^C}^{n+1} = T_{I^C}^n - v_x \Delta{t}\frac{T_{I^E}^n - T_{I^W}^n}{2\Delta{x}}.
\end{equation}$

The right-hand side can be simplified using the so-called *Courant number*:

$\begin{equation}
\alpha = \frac{v_x\Delta{t}}{\Delta{x}},
\end{equation}$

which represents the number of grid points traversed in a single time step. 

Unfortunately, this scheme is unconditionally unstable for the advection equation, as shown by a *Von Neumann* or *Hirt's stability analysis*. The central difference at $I^C$ causes amplification of the variable (here, temperature) at each subsequent time step. Hence, the solution continually grows and is unstable.

## Lax-Friedrichs method

One way to suppress the instability of the FTCS scheme is the *Lax-Friedrichs* method. This replaces the term $T_{I^C}^{n}$ with its spatial average at the same time level, resulting in:

$\begin{equation}
\frac{T_{I^C}^{n+1}-\left(T_{I^E}^{n}+T_{I^W}^{n}\right)/2}{\Delta{t}}=-v_x\frac{T_{I^E}^{n}-T_{I^W}^{n}}{2\Delta{x}}.
\end{equation}$

Rearanging gives:

$\begin{equation}
T_{I^C}^{n+1} = \frac{1}{2}\left(T_{I^E}^{n}+T_{I^W}^{n}\right)-
\frac{v_x \Delta{t}}{2\Delta{x}} \left(T_{I^E}^{n}-T_{I^W}^{n}\right).
\end{equation}$

This method is stable for $\alpha < 1$ but introduces significant numerical diffusion.

## Upwind

Another approach is to consider only upstream information. The *upwind* scheme uses one-sided finite differences, always taken in the upstream direction. This results in a scheme that is first-order accurate in both space and time. The discretized advection equation becomes:

$\begin{equation}
\frac{T_{I^C}^{n+1}-T_{I^C}^n}{\Delta{t}} = -v_{x,I^C}
\begin{cases}
\frac{T_{I^C}^{n}-T_{I^W}^{n}}{\Delta{x}} &\text{if } v_{x,I^C} \gt 0\\
\frac{T_{I^E}^{n}-T_{I^C}^{n}}{\Delta{x}}&\text{if } v_{x,I^C} \lt 0 
\end{cases}.
\end{equation}$

The scheme is stable if the CFL-criterion is satisfied ($\alpha \le 1$), but numerical diffusion remains, which depends on the grid size. A Taylor series expansion shows that, in 1D with constant velocity, the scheme becomes non-diffusive if the time step exactly satisfies the CFL-criterion. The method becomes unstable if this criterion is violated.

## Staggered Leapfrog

All previously discussed explicit schemes are only first-order accurate in time and second-order in space (except upwind, which is first-order in both). To match the temporal and spatial accuracy without choosing a very small time step, one may use the *staggered leapfrog* scheme:

$\begin{equation}
\frac{T_{I^C}^{n+1}-T_{I^C}^{n-1}}{2\Delta{t}}=-v_x\frac{T_{I^E}^{n}-T_{I^W}^{n}}{2\Delta{x}}.
\end{equation}$

This method avoids numerical diffusion, but becomes increasingly unstable when strong gradients in the advected field are present.

## Semi-Lagrangian

The methods discussed above each have drawbacks. The *semi-Lagrangian* method addresses several of them: it is stable, does not suffer from numerical diffusion, and is not constrained by the CFL criterion. It is related to tracer-based advection schemes and solves ODEs rather than using traditional finite differences. While not inherently conservative and subject to minor interpolation errors, it offers promising accuracy and efficiency.

The central idea is to trace an advected particle backward in time to its origin and interpolate the corresponding value from the Eulerian grid.

In 1D, assuming constant velocity in time and space, the procedure is:

**1. Calculate the initial position** 

The initial position $X_{i}$ of a particle landing on the Eulerian grid point $x_{I^C}$ at time step ${n+1}$ is:

$\begin{equation}
X_{i}=x_{I^C}-\Delta{t}\cdot v_{x,I^C}^{n+1},
\end{equation}$

where $x_{I^C}$ is the coordinate of the Eulerian grid point $I^C$, $\Delta{t}$ is the time step, $v_x$ the velocity in $x$-direction, and $n+1$ is the time at the new time step. 

**2. Interpolate the temperature**

Interpolate the temperature at time step $n$ from the surrounding Eulerian grid points onto the position $X_{i}$, e.g., using `cubic_spline_interpolation()`.

**3. Update the temperature field**

Assuming the temperature at the grid point $I^C$ at time step $n+1$ equals the interpolated value at $X_{i}$ at time step $n$ results in:

$\begin{equation}
T_{I^C}^{n+1} = T_{X_{I^C}}^n
\end{equation}$

<!-- Stop -->

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


# Advection Equation (1D)

In one spatial dimension, the advection equation for temperature (assuming incompressible flow) is given by

$\begin{equation}
\frac{\partial{T}}{\partial{t}} = -v_x \left(\frac{\partial{T}}{\partial{x}}\right),
\end{equation}$

where $T$ denotes temperature [K], $t$ is time [s], and $v_x$ is the velocity in the $x$-direction. In the following 1D derivations, a spatially and temporally constant velocity $v_x$ is assumed unless stated otherwise.

# Discretization Schemes

The global indexing of the central reference point $I$ follows the convention introduced in the [general solution section](./GESolution.md). The adjacent indices are defined as

$\begin{equation}
\begin{split}
I^\textrm{W} &= I^\textrm{C} - 1, \\
I^\textrm{C} &= I, \\
I^\textrm{E} &= I^\textrm{C} + 1,
\end{split}
\end{equation}$

where $I$ denotes the equation number corresponding to the local grid index $i$. The superscripts $C$, $W$, and $E$ indicate the central, western, and eastern grid points of the three-point stencil used below.

## Forward in Time, Centered in Space (FTCS)

Approximating the temporal derivative with a forward difference and the spatial derivative with a centered difference yields

$\begin{equation}
\frac{T_{I^\textrm{C}}^{n+1}-T_{I^\textrm{C}}^{n}}{\Delta{t}} = -v_x\left(\frac{T_{I^\textrm{E}}^{n}-T_{I^\textrm{W}}^{n}}{2\Delta{x}}\right),
\end{equation}$

where $\Delta{t}$ and $\Delta{x}$ denote the time step and grid spacing, respectively, and $n$ is the temporal index. The scheme is first-order accurate in time and second-order accurate in space.

Rearranging gives

$\begin{equation}
T_{I^\textrm{C}}^{n+1} = T_{I^\textrm{C}}^n - v_x \Delta{t}\frac{T_{I^\textrm{E}}^n - T_{I^\textrm{W}}^n}{2\Delta{x}}.
\end{equation}$

Introducing the Courant number

$\begin{equation}
\alpha = \frac{v_x\Delta{t}}{\Delta{x}},
\end{equation}$

which represents the fraction of a grid cell traversed during one time step. For explicit finite-difference schemes, stability typically requires satisfaction of the Courant–Friedrichs–Lewy (CFL) condition $|\alpha| \le 1$, ensuring that information does not propagate more than one grid cell per time step.

Unfortunately, the FTCS scheme is unconditionally unstable for the advection equation, as shown by a von Neumann stability analysis. The centered spatial discretization at $I^\textrm{C}$ leads to amplification of perturbations, resulting in an unstable solution.

## Lax-Friedrichs Method

The Lax–Friedrichs method stabilizes the FTCS scheme by replacing $T_{I^\textrm{C}}^{n}$ with its spatial average at the same time level:

$\begin{equation}
\frac{T_{I^\textrm{C}}^{n+1}-\left(T_{I^\textrm{E}}^{n}+T_{I^\textrm{W}}^{n}\right)/2}{\Delta{t}}=-v_x\frac{T_{I^\textrm{E}}^{n}-T_{I^\textrm{W}}^{n}}{2\Delta{x}}.
\end{equation}$

Rearranging gives

$\begin{equation}
T_{I^\textrm{C}}^{n+1} = \frac{1}{2}\left(T_{I^\textrm{E}}^{n}+T_{I^\textrm{W}}^{n}\right)-
\frac{v_x \Delta{t}}{2\Delta{x}} \left(T_{I^\textrm{E}}^{n}-T_{I^\textrm{W}}^{n}\right).
\end{equation}$

The scheme is stable for $\alpha \le 1$ but introduces significant numerical diffusion.

## Upwind Scheme

The upwind scheme accounts for the direction of information propagation by using one-sided spatial differences. The discretized equation reads

$\begin{equation}
\frac{T_{I^\textrm{C}}^{n+1}-T_{I^\textrm{C}}^n}{\Delta{t}} = -v_{x}
\begin{cases}
\frac{T_{I^\textrm{C}}^{n}-T_{I^\textrm{W}}^{n}}{\Delta{x}} &\text{if } v_{x} \gt 0\\
\frac{T_{I^\textrm{E}}^{n}-T_{I^\textrm{C}}^{n}}{\Delta{x}}&\text{if } v_{x} \lt 0 
\end{cases}.
\end{equation}$

The scheme is first-order accurate in both time and space and is conditionally stable, requiring satisfaction of the CFL condition ($|\alpha| \le 1$). However, it introduces numerical diffusion proportional to the grid spacing. For constant velocity in 1D, the scheme becomes non-diffusive if the time step exactly satisfies the CFL condition.

---

## Staggered Leapfrog

To achieve second-order accuracy in both space and time, the staggered leapfrog scheme can be used:

$\begin{equation}
\frac{T_{I^\textrm{C}}^{n+1}-T_{I^\textrm{C}}^{n-1}}{2\Delta{t}}=-v_x\frac{T_{I^\textrm{E}}^{n}-T_{I^\textrm{W}}^{n}}{2\Delta{x}}.
\end{equation}$

This method is non-diffusive but can exhibit dispersive oscillations, particularly near sharp gradients.

## Semi-Lagrangian Method

The semi-Lagrangian method combines Eulerian and Lagrangian concepts. It is unconditionally stable with respect to the CFL condition and significantly reduces numerical diffusion, though it is not strictly conservative and depends on interpolation accuracy.

The central idea is to trace a particle backward in time to determine its departure point and interpolate the corresponding value from the grid.

Assuming constant velocity in 1D:

**1. Calculate the initial position** 

The initial position $X_{i}$ of a particle landing on the Eulerian grid point $x_{I^\textrm{C}}$ at time step ${n+1}$ is:

$\begin{equation}
X_{i}=x_{I^\textrm{C}}-\Delta{t}\cdot v_{x},
\end{equation}$

where $x_{I^\textrm{C}}$ is the coordinate of the Eulerian grid point $I^\textrm{C}$. 

**2. Interpolate the temperature**

Interpolate $T^n$ from the surrounding grid nodes onto $X_i$ (e.g., using cubic spline interpolation).

**3. Update the temperature**

Assuming the temperature at the grid point $I^\textrm{C}$ at time step $n+1$ equals the interpolated value at $X_{i}$ at time step $n$ results in:

$\begin{equation}
T_{I^\textrm{C}}^{n+1} = T^n(X_{i})
\end{equation}$

## Passive Tracers

Passive tracers represent a fully Lagrangian approach. Tracers with initial positions $x_p(t=0)$ and associated properties are transported by solving the particle trajectory equation

$\begin{equation}
\frac{dx_p}{dt}=v_x \left(x_p,t\right).
\end{equation}$

**Forward Euler**

The flow path ODE is approximated as:

$\begin{equation}
\frac{x_p^{n+1}-x_p^n}{\Delta{t}} = v_x. 
\end{equation}$

Solving for the next position:

$\begin{equation}
x_p^{n+1} = x_p^n + \Delta{t}\cdot v_x. 
\end{equation}$

This method is simple but inaccurate for large time steps.
    
**Fourth-Order Runge-Kutta**

In 1D, the next position is:

$\begin{equation}
x_p^{n+1} = x_p^n + \frac{1}{6}k_1 + \frac{1}{3}k_2 + \frac{1}{3}k_3 + \frac{1}{6}k_4,
\end{equation}$

with

$\begin{equation}
\begin{split}
k_1 & = \Delta{t} \cdot v_x \\
k_2 & = \Delta{t} \cdot v_x \\
k_3 & = \Delta{t} \cdot v_x \\
k_4 & = \Delta{t} \cdot v_x \\
\end{split}
\end{equation}$

Tracer properties are subsequently interpolated back to the Eulerian grid as required, e.g., using linear interpolation.

> Note: For constant velocity, the fourth-order Runge–Kutta scheme reduces to the Forward Euler update.

While highly flexible, tracer methods require careful treatment. Interpolation between grid and tracer data may introduce smoothing, and tracer clustering or depletion can lead to reduced accuracy. Currently, adaptive tracer correction techniques are not implemented in `GeoModBox.jl`.
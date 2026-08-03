# Advection Equation (1D)

In one spatial dimension, assuming a spatially and temporally constant advection
velocity, the governing advection equation is given by

$\begin{equation}
\frac{\partial T}{\partial t} = -v_x \left(\frac{\partial T}{\partial x}\right),
\end{equation}$

where $T$ denotes temperature [K], $t$ is time [s], and $v_x$ is the velocity
in the $x$-direction. Throughout the following derivations, the advection
velocity is assumed to be spatially and temporally constant unless stated
otherwise.

# Discretization Schemes

The global indexing of the central reference point $I$ follows the convention
introduced in the [general solution section](./GESolution.md). The adjacent
indices are defined as

$\begin{equation}
\begin{split}
I^\textrm{W} &= I^\textrm{C} - 1, \\
I^\textrm{C} &= I, \\
I^\textrm{E} &= I^\textrm{C} + 1,
\end{split}
\end{equation}$

where $I$ denotes the equation number corresponding to the local grid index
$i$. The superscripts $C$, $W$, and $E$ indicate the central, western, and
eastern grid points of the three-point stencil used throughout this chapter.

## Forward in Time, Centered in Space (FTCS)

Approximating the temporal derivative with a forward difference and the
spatial derivative with a centered difference yields

$\begin{equation}
\frac{T_{I^\textrm{C}}^{n+1}-T_{I^\textrm{C}}^{n}}{\Delta t} =
-v_x \left(
\frac{T_{I^\textrm{E}}^{n}-T_{I^\textrm{W}}^{n}}
{2\Delta x}
\right),
\end{equation}$

where $\Delta t$ and $\Delta x$ denote the time step and grid spacing,
respectively, and $n$ is the temporal index. The scheme is first-order
accurate in time and second-order accurate in space.

Rearranging gives

$\begin{equation}
T_{I^\textrm{C}}^{n+1} =
T_{I^\textrm{C}}^n -
v_x \Delta t
\frac{T_{I^\textrm{E}}^n-T_{I^\textrm{W}}^n}{2\Delta x}.
\end{equation}$

Introducing the Courant number

$\begin{equation}
\alpha=\frac{v_x\Delta t}{\Delta x},
\end{equation}$

which represents the fraction of a grid cell traversed during one time step.
Although explicit finite-difference schemes generally require $|\alpha|\le1$,
the FTCS scheme is unconditionally unstable for the advection equation.

## Lax-Friedrichs Method

The Lax-Friedrichs method stabilizes the FTCS scheme by replacing
$T_{I^\textrm{C}}^{n}$ with its spatial average:

$\begin{equation}
\frac{T_{I^\textrm{C}}^{n+1}-\left(T_{I^\textrm{E}}^{n}+T_{I^\textrm{W}}^{n}\right)/2}{\Delta t} =
-v_x\frac{T_{I^\textrm{E}}^{n}-T_{I^\textrm{W}}^{n}}{2\Delta x}.
\end{equation}$

Rearranging gives

$\begin{equation}
T_{I^\textrm{C}}^{n+1} =
\frac{1}{2}\left(T_{I^\textrm{E}}^{n}+T_{I^\textrm{W}}^{n}\right) -
\frac{v_x\Delta t}{2\Delta x}
\left(T_{I^\textrm{E}}^{n}-T_{I^\textrm{W}}^{n}\right).
\end{equation}$

The scheme is conditionally stable for $|\alpha|\le1$, but introduces
significant numerical diffusion owing to the additional averaging step.

## Upwind Scheme

The upwind scheme accounts for the direction of information propagation by
using one-sided spatial differences:

$\begin{equation}
\frac{T_{I^\textrm{C}}^{n+1}-T_{I^\textrm{C}}^n}{\Delta t} =
-v_x
\begin{cases}
\frac{T_{I^\textrm{C}}^{n}-T_{I^\textrm{W}}^{n}}{\Delta x} &\text{if } v_x>0\\
\frac{T_{I^\textrm{E}}^{n}-T_{I^\textrm{C}}^{n}}{\Delta x} &\text{if } v_x<0
\end{cases}.
\end{equation}$

The scheme is first-order accurate in both space and time and is conditionally
stable for $|\alpha|\le1$. Numerical diffusion depends on the Courant number
and decreases with increasing spatial and temporal resolution. For a constant
velocity field in one dimension, the scheme becomes non-diffusive when the
Courant number is exactly equal to one.

## Staggered Leapfrog

To achieve second-order accuracy in both space and time, the staggered
leapfrog scheme can be used:

$\begin{equation}
\frac{T_{I^\textrm{C}}^{n+1}-T_{I^\textrm{C}}^{n-1}}{2\Delta t} =
-v_x
\frac{T_{I^\textrm{E}}^{n}-T_{I^\textrm{W}}^{n}}{2\Delta x}.
\end{equation}$

The scheme exhibits very little numerical diffusion but may develop dispersive
oscillations near sharp gradients. Because the leapfrog method requires the
solution at two previous time levels, the first time step must be initialized
using a separate scheme.

## Semi-Lagrangian Method

The semi-Lagrangian method combines Eulerian and Lagrangian concepts. The
method is not restricted by the conventional explicit CFL stability condition.
However, the time step still affects the accuracy of the calculated
trajectories and the interpolation of the transported property.

**1. Determine the departure point**

$\begin{equation}
X_i=x_{I^\textrm{C}}-\Delta t\,v_x,
\end{equation}$

**2. Interpolate the temperature**

The temperature field at time level $n$ is interpolated onto the departure
point $X_i$. In the current implementation of `GeoModBox.jl`, linear spline
interpolation is used.

**3. Update the temperature**

$\begin{equation}
T_{I^\textrm{C}}^{n+1}=T^n(X_i).
\end{equation}$

## Passive Tracers

Passive tracers provide a fully Lagrangian description of material transport.

$\begin{equation}
\frac{dx_p}{dt}=v_x(x_p,t).
\end{equation}$

**Forward Euler**

$\begin{equation}
x_p^{n+1}=x_p^n+\Delta t\,v_x.
\end{equation}$

**Fourth-Order Runge-Kutta**

For a spatially and temporally constant velocity field, all four Runge-Kutta
stages evaluate the same velocity,

$\begin{equation}
k_1=k_2=k_3=k_4=\Delta t\,v_x,
\end{equation}$

such that the fourth-order Runge-Kutta method reduces exactly to

$\begin{equation}
x_p^{n+1}=x_p^n+\Delta t\,v_x.
\end{equation}$

Although highly flexible, tracer methods require careful treatment.
Interpolation between the Eulerian grid and the tracer distribution may
introduce smoothing errors, while tracer clustering or depletion can reduce
the local interpolation accuracy. Adaptive tracer insertion or redistribution
techniques are currently not implemented in `GeoModBox.jl`.

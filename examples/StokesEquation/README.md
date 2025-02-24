# Momentum Conservation Equation

<!-- -->

------------------
------------------

## Examples
<!--
### Channel Flow
&emsp; Assuming the horizontal pressure gradient is constant and flow within a channel is only driven by the pressure and/or by a constant horizonal velocity at the surface (or at the bottom, or both), the stokes equation describes the horizontal flow velocity within the channel and simplifies to: 

$$
\frac{\partial P}{\partial x} = \frac{\partial \tau_{x,z}}{\partial z}
$$,

where *P* is the pressure and $\tau_{x,z}$ is the deviatoric shear stress, which is defined as: 

$\tau_{x,z} = 2 \eta \dot{\varepsilon}_{x,z} = \eta \frac{\partial v_x}{\partial z}$,

where $\eta$ is the dynamic viscosity and $\dot{\varepsilon}_{x,z}$ is the deviatoric shear strain-rate, which is defined as: 

$\dot{\varepsilon}_{x,z} = \frac{1}{2}(\frac{\partial v_x}{\partial z} + \frac{\partial v_z}{\partial x})$.

For the given setup I can assume that the vertical velocity is zero and thus equation (3) simplifies to the last expression of equation (2).

&emsp; This directory contains a script to calculate the horizontal velocity for a two-dimensional Couette(-Poiseuille) channel flow with constant and logarithmically, with depth varying viscosity and to compare the numerical solution with its analytical solution. The depth-dependent viscosity is defined as: 

$\eta = \eta_0 exp(log(m) \frac{H-z}{H})$,

where *m* is the viscosity ratio of $\frac{\eta_1}{\eta_0}$, $\eta_0$ and $\eta_1$ are the bottom and surface viscosities, respectively, *H* is the model height, and *z* the depth. 

&emsp;Considering the definition of the viscosity as given in equation (4), one can derive an analytical solution of the horizontal velocity from the 1-D stokes equation in *x*-direction by twice integrating equation (1). The analytical solution with depth depends on the viscosity ratio *m*, the horizontal pressure gradient $\frac{\partial P}{\partial x}$, and the shear velocity at the surface $v_{x,0}$. For an upward pointing coordinate system (*z* positive) the analytical solution is given as: 

$v_{x,ana} =-\frac{1}{2 \eta_0} \frac{\partial P}{x} (Hz - z^2) + v_{x,0}\frac{z}{H}$,&emsp;&emsp; if $m = 1$, and &emsp;&emsp;&emsp; (5)

$v_{x,ana} = -\frac{\partial P}{\partial x} \frac{H}{\eta_0 log(m)} (\frac{m^{-\frac{z}{H}}}{m-1}(z(m-1)+H) - \frac{H}{m-1})-m^{-\frac{z}{H}} m \frac{v_{x,0}}{m-1} + \frac{v_{x,0}m}{m-1}$, &emsp;&emsp; if $m \neq 0$.&emsp;&emsp;&emsp; (6)

&emsp;The numerical solution is calculated using fixed boundary velocities, which are defined by the analytical solution of the horizontal velocity as defined in equations (5) and (6) and I simply flip the analytical solution so that it fits to the downward point coordinate system I use in the code. -->

---------------------------------------------------------------

<!-- Channel flow example in quasi 2-D:
    - prescribe velocity boundary conditions by the analytical solution of a channel flow 
-->
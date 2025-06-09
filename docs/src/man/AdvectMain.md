# Advection Equation

In case the material is not moving, one can solve the energy equation only for the diffusive part (e.g., an intrusion problem or a non-deforming lithosphere). Generally, however, the material is moving and certain properties need to be advected with the flow (e.g., the temperature, density, composition etc.). 
   
> **Note:** The energy equation can be solved simultaneously with the diffusive and convective part using different discretization methods (interestingly, *forward in time and centered in space (FTCS)* is stable with some numerical diffusion, which is always unstable for pure advection). 

For the sake of simplicity and a more conveniant way to teach both mechanisms, so far, the operator-splitting method is preferred within the ```GeoModBox.jl```. Within the operator-splittin method, first the convective part of the temperature conservation equation is solved, followed by the conductive part. The conducitve part can be solved by different discretization methods as described in the [diffusion equation documentation](./DiffMain.md). For the convective part, the ```GeoModBox.jl``` focuses on the

- upwind,
- lax, 
- staggered-leaped frog,
- semi-lagrangian, and
- tracer method. 

In general, advection describes the transport of a property from one point to another, where one can assume different reference frames for the given point of interest. If one assumes a not moving reference frame (that is an *Eulerian* grid), the change in temperature at a certain point can be described by (i.e. the *eulerian advective transport equation*)

$\begin{equation}
\frac{\partial T}{\partial t} = - \overrightarrow{v} \cdot \nabla{T},
\end{equation}$

or in a Lagrangian reference frame (along a moving point; i.e., the *substantive* derivative) as

$\begin{equation}
\frac{DT}{Dt},
\end{equation}$

where both are related by

$\begin{equation}
\frac{DT}{Dt} = \frac{\partial T}{\partial t} + \overrightarrow{v} \cdot \nabla{T}.
\end{equation}$

For a Lagrangian reference point, advection is given by a simple *ordinary differential equation* particle advection scheme, where changes in its coordinates are related with the material velocities as

$\begin{equation}
\frac{Dx_i}{Dt} = v_i,
\end{equation}$

where $i$ is the coordinate index and $x_i$ is a spatial coordinate. 

## Discretization Schemes

As simple as the advection equation seems to be, it is rather difficult to properly solve advection without some kind of numerical diffusion or inaccuracies due to interpolation of properties between the tracers and the regular grid. So care needs to be taken.  

The particle advection is used for the *tracer/marker in cell* method (either passive or active) and can be solved using different numerical methods, e.g., Euler integration or Runge-Kutta. The Eulerian form of the advection equation can also be solved in different ways.

<!-- The ```GeoModBox.jl``` focus on *four* different methods to advect material. Different advection methods are used within the individual benchmarks and all can be tested in the [*Rigid Body Rotation*](https://github.com/LukasFuchs/FDCSGm/tree/main/Benchmark/RigidBodyRotation) benchmark.  -->

**Advection examples:**

- [2-D advection with constant velocity field](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/AdvectionEquation/2D_Advection.jl)  
- [Resolution test of 2-D advection](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/AdvectionEquation/2D_Advection_ResolutionTest.jl)

See the [examples documentation](./Examples.md) for further details.

**Advection exercises include:**

- [1-D Gaussian or block anomaly advection](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/06_1D_Advection.ipynb)  
- [2-D coupled advection-diffusion](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/07_2D_Energy_Equation.ipynb)

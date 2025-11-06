# Specifications

`GeoModBox.jl` includes several [routines](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/InitialCondition/2Dini.jl) or [structures](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/Structures.jl) to define certain parameters or initialize specific anomalies. The initial conditions can be specified for properties defined on their corresponding grid (i.e., temperature, velocity, or phase) or for tracers. 

## Initial Temperature

```julia
IniTemperature!(type,M,NC,D,x,y;Tb=600.0,Ta=1200.0,σ=0.1)
```

Currently, there are five different initial temperature conditions available: 

1. A circular anomaly with a constant background (`circle`)
2. A Gaussian anomaly (`gaussian`)
3. A rectangular shaped anomaly with a constant background (`block`)
4. A linear increasing temperature with depth (`linear`)
5. A linear increasing temperature with depth including an elliptical anomaly (`lineara`)

The input parameters are: 

- type - Parameter defining the type (see above)
- M - Structure or tuple containing the geometry
- NC - Structure or tuple containing the centroids parameter 
- D - Structure or tuple containing the field arrays
- x - Structure or tuple containing the x-coordinates
- y - Structure or tuple containing the y-coordinates

Certain default values can be modified as well: 

- Tb - Scalar value for the background (or top) temperature
- Ta - Scalar value for the maximum (anomaly or bottom) temperature
- σ - Width of the Gaussian temperature anomaly

The temperature is initialized on the extended centroid grid. The corresponding field without *ghost nodes* is updated accordingly.

The function is called, for example, like [here](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/MixedHeatedConvection/BottomHeated.jl): 

```Julia 
IniTemperature!(Ini.T,M,NC,D,x,y;Tb=P.Tbot,Ta=P.Ttop)
``` 

## Initial Velocity

```julia
IniVelocity!(type,D,VBC,NV,Δ,M,x,y;ε=1e-15)
```

The following velocity configurations are currently supported: 

1. A rigid-body rotation (`RigidBody`)
2. A shear cell (`ShearCell`)
3. Simple Shear (`SimpleShear`) 
4. Pure Shear (`PureShear`)

The input parameters are: 

- type - Parameter defining the type (see above)
- D - Structure or tuple containing the field arrays
- VBC - Structure or tuple containing the velocity boundary conditions
- NV - Structure or tuple containing the vertices parameter 
- Δ - Structure or tuple containing the grid resolution
- M - Structure or tuple containing the geometry
- x - Structure or tuple containing the x-coordinates
- y - tructure or tuple containing the y-coordinates

Certain default values can be modified as well: 

- ε - Background strain rate for pure shear or simple shear

The function is called, for example, like [here](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/AdvectionEquation/2D_Advection.jl): 

```Julia
IniVelocity!(Ini.V,D,VBC,NC,NV,Δ,M,x,y)
``` 

## Initial Phase

```julia
IniPhase!(type,D,M,x,y,NC;phase=0)
```

Currently, only one initial phase configuration is available: 

1. A rectangular shaped anomaly (`block`)

The input parameters are: 

- type - Parameter defining the type (see above)
- D - Structure or tuple containing the field arrays
- M - Structure or tuple containing the geometry
- x - Structure or tuple containing the x-coordinates
- y - Structure or tuple containing the y-coordinates
- NC - Structure or tuple containing the centroids parameter

Certain default values can be modified as well: 

- phase - Vector containing the phase ID numbers (e.g. phase=[0,1])

The phase is initialized on the extended centroid grid. The corresponding field without *ghost nodes* is updated accordingly.

The function is called, for example, like [here](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/StokesEquation/2D/FallingBlockBenchmark.jl): 

```Julia
IniPhase!(Ini.p,D,M,x,y,NC;phase)
```

The density on the extended, centroid grid is then, for example, updated via the phase ID: 

```julia
for i in eachindex(phase)
    D.ρ[D.p.==phase[i]] .= ρ[i]
end
```

# Tracer Method

Because tracer advection can be performed in parallel, additional parameters must be defined. `GeoModBox.jl` provides functionality to initialize certain tracer positions. Additional initial configuration methods are encouraged and can be integrated. The properties can either be interpolated from the centroids to the tracers or initialized first on the tracers followed by an interpolation on the grid. 

The following steps are required to use tracers: 

**1. Tracer initialization**

To initialize the tracers, one needs to define what property should be advected, the number of markers per cell, the wanted noise, the initial shape, and the phase ID. The tracer properties can be interpolated to the centroids or vertices using a bilinear interpoltion scheme. For the centroids, the extended centroid field must be used. The remaining parameters are the general `tuples` or `structures` used in `GeoModBox.jl`.

Following the definition of the required parameters for the tracer advection, the initial tracer position can be defined via the function 

```julia
 IniTracer2D(Aparam,nmx,nmy,Δ,M,NC,noise,ini,phase;λ=1.0e3,δA=5e2/15,ellA=100.0,ellB=100.0,α=0.0)
``` 

The function initializes the position, phase, and memory of the tracers. As initial tracer phase distribution one can choose: 
-   `ini=:block`        : a rectangular block
-   `ini=:RTI`          : a cosine perturbation with wavelength λ and amplitude δA 
    `ini:=Inclusion`    : a viscous circular or elliptical inclusion

The input parameters are: 

    - Aparam - defines if temperature (`thermal`) or phase (`phase`) is advected
    - nmx - number of horizontal tracers per cell
    - nmy - number of vertical tracers per cell
    - Δ - Structure or tuple containing the grid resolution
    - M - Structure or tuple containing the geometry
    - NC - Structure or tuple containing centroids parameter
    - noise - add noise; 1 - yes, 0 - no
    - ini - Initial phase distribution (`block`)
    - phase - Vector with phase IDs, (e.g. [0,1])

Certain default values can be modified as well: 

    - λ - Wavelength [m] for a cosine perturbation, e.g. for the RTI
    - δA - Amplitude [m] of the perturbation
    - ellA - Major half axis [m] of the elliptical inclusion
    - ellB - Minor half axis [m] of the elliptical inclusion
    - α - Rotation angle [°] of the elliptical inclusion  

The output is a tuple including the x-, and y-coordinates of the tracers and the phase ID as one dimensional arrays. If temperature is advected, temperature is also stored for each marker. 

To advect the temperature, the initialization is called, for example, like [here](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/AdvectionEquation/2D_Advection.jl): 

```julia
# Tracer Initialization ---
nmx,nmy     =   3,3
noise       =   1
nmark       =   nmx*nmy*NC.x*NC.y
Aparam      =   :thermal
MPC         =   (
    c       =   zeros(Float64,(NC.x,NC.y)),
    v       =   zeros(Float64,(NV.x,NV.y)),
    th      =   zeros(Float64,(nthreads(),NC.x,NC.y)),
    thv     =   zeros(Float64,(nthreads(),NV.x,NV.y)),
)
MAVG        = (
    PC_th   =   [similar(D.wte) for _ = 1:nthreads()],  # per thread
    PV_th   =   [similar(D.wtv) for _ = 1:nthreads()],   # per thread
    wte_th  =   [similar(D.wte) for _ = 1:nthreads()],  # per thread
    wtv_th  =   [similar(D.wtv) for _ = 1:nthreads()],  # per thread
)
# MPC     =   merge(MPC,MPC1)
Ma      =   IniTracer2D(Aparam,nmx,nmy,Δ,M,NC,noise,0,0)
# RK4 weights ---
rkw     =   1.0/6.0*[1.0 2.0 2.0 1.0]   # for averaging
rkv     =   1.0ftz/2.0*[1.0 1.0 2.0 2.0]   # for time stepping
# Interpolate on centroids ---
@threads for k = 1:nmark
    Ma.T[k] =   FromCtoM(D.T_ex, k, Ma, x, y, Δ, NC)
end
# Count marker per cell ---
CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,NV,1)
```

To advect the phase, the initialization is called, for example, like [here](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/StokesEquation/2D/FallingBlockVarEta_DC.jl): 

```julia 
# Tracer Advection ================================================== #
nmx,nmy     =   3,3
noise       =   0
nmark       =   nmx*nmy*NC.x*NC.y
Aparam      =   :phase
MPC         =   (
        c       =   zeros(Float64,(NC.x,NC.y)),
        v       =   zeros(Float64,(NV.x,NV.y)),
        th      =   zeros(Float64,(nthreads(),NC.x,NC.y)),
        thv     =   zeros(Float64,(nthreads(),NV.x,NV.y)),
)
MAVG        = (
        PC_th   =   [similar(D.wte) for _ = 1:nthreads()],  # per thread
        PV_th   =   [similar(D.ηv) for _ = 1:nthreads()],   # per thread
        wte_th  =   [similar(D.wte) for _ = 1:nthreads()],  # per thread
        wtv_th  =   [similar(D.wtv) for _ = 1:nthreads()],  # per thread
)
# MPC     =   merge(MPC,MPC1)
Ma      =   IniTracer2D(Aparam,nmx,nmy,Δ,M,NC,noise,Ini.p,phase)
# RK4 weights ---
rkw     =   1.0/6.0*[1.0 2.0 2.0 1.0]   # for averaging
rkv     =   1.0/2.0*[1.0 1.0 2.0 2.0]   # for time stepping
# Count marker per cell ---
CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,NV,1)
# Interpolate from markers to cell ---
Markers2Cells(Ma,nmark,MAVG.PC_th,D.ρ_ex,MAVG.wte_th,D.wte,x,y,Δ,Aparam,ρ)
D.ρ     .=  D.ρ_ex[2:end-1,2:end-1]  
Markers2Cells(Ma,nmark,MAVG.PC_th,D.p_ex,MAVG.wte_th,D.wte,x,y,Δ,Aparam,phase)
D.p     .=  D.p_ex[2:end-1,2:end-1]
Markers2Cells(Ma,nmark,MAVG.PC_th,D.η_ex,MAVG.wte_th,D.wte,x,y,Δ,Aparam,η)
D.ηc    .=  D.η_ex[2:end-1,2:end-1]
Markers2Vertices(Ma,nmark,MAVG.PV_th,D.ηv,MAVG.wtv_th,D.wtv,x,y,Δ,Aparam,η)
```

To interpolate a property, like the temperature, from the centroids to the tracers one can use the function `FromCtoM()`. 

To interpolate a property from the tracers back to the centroids or vertices, like the viscosity, one can use the functions `Markers2Cells()` or `Markers2Vertices()`. Per default, an arithmetic averaging scheme is used to inerpolate the property, however, geometric and harmonic means are also available. For more information regarding the functions, please see the help function. 

**2. Tracer advection** 

The tracers are advected using Runge-Kutta 4th order. This is conducted using the function

```julia
AdvectTracer2D(Ma,nmark,D,x,y,dt,Δ,NC,rkw,rkv,style)
```

The input parameters are: 

- Ma - Structure containing the tracer information
- nmark - Total number of tracers
- D - Structure or tuple containing the field arrays
- x - Structure or tuple containing the x-coordinates
- y - Structure or tuple containing the y-coordinates
- dt - Time step
- Δ - Structure or tuple containing the grid resolution
- NC - Structure or tuple containing the centroids parameter
- rkw - Runge-Kutta weights for averaging
- rkv - Runge-Kutta weights for time stepping
- style - Defines which velocity is used for the advection

For more details please refer to the [source code](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/Tracers/2Dsolvers.jl).

>**Note:** Currently, temperature is not intended to be advected via tracers, as this would require the update of the tracer temperature via incremental changes rather than the absolute value. Within the [2-D advection example](./examples/Advection2D.md) temperature advection is only used assuming non-diffusive process. Thus, no update of the tracer temperature is required! 

The advection of temperature and the update of the temperature field on the centroids is called, for example, like [here](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/AdvectionEquation/2D_Advection.jl): 

```julia 
# Advect tracers ---
AdvectTracer2D(Ma,nmark,D,x,y,T.Δ[1],Δ,NC,rkw,rkv,1)
# CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,i)
CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,NV,i)
     
# Interpolate temperature from tracers to grid ---
Markers2Cells(Ma,nmark,MAVG.PC_th,D.T_ex,MAVG.wte_th,D.wte,x,y,Δ,Aparam,0)           
D.T     .=  D.T_ex[2:end-1,2:end-1]
```

The advection of the phase and the update of the corresponding grid parameters is called, for example, like [here](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples/StokesEquation/2D/FallingBlockVarEta_DC.jl): 

```julia
# Advection ===
# Advect tracers ---
@printf("Running on %d thread(s)\n", nthreads())  
AdvectTracer2D(Ma,nmark,D,x,y,T.Δ[1],Δ,NC,rkw,rkv,1)
CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,NV,it)

# Interpolate phase from tracers to grid ---
Markers2Cells(Ma,nmark,MAVG.PC_th,D.ρ_ex,MAVG.wte_th,D.wte,x,y,Δ,Aparam,ρ)
D.ρ     .=   D.ρ_ex[2:end-1,2:end-1]  
Markers2Cells(Ma,nmark,MAVG.PC_th,D.p_ex,MAVG.wte_th,D.wte,x,y,Δ,Aparam,phase)
D.p     .=  D.p_ex[2:end-1,2:end-1]
Markers2Cells(Ma,nmark,MAVG.PC_th,D.η_ex,MAVG.wte_th,D.wte,x,y,Δ,Aparam,η)
D.ηc    .=   D.η_ex[2:end-1,2:end-1]
Markers2Vertices(Ma,nmark,MAVG.PV_th,D.ηv,MAVG.wtv_th,D.wtv,x,y,Δ,Aparam,η)
```

>**Note:** The tracer distribution and interpolation of tracer properties to the centroids or vertices is a very helpful feature to initialize different, more complex model setups. This will be part of future implementations. 

# Mutable Structures

Implementation of scaled thermal convection models in `GeoModBox.jl` revealed that immutable `NamedTuples` can be numerically limiting. Some initially defined parameters, like the model height, need to be modified in order to scale them. Thus, mutable structures (`mutable struct`) have been included. Luckily, the functions do not distinguish between a `NamedTuple` and a `mutable struct` in `Julia`. Thus, no additional modification in the functions are required. 

Mutable structures are designed so that all necessary parameters are initialized with default values. To ensure that the default values are picked if no input parameter is used calling the `mutable struct` the option `@kwdef` needs to be added when initializing the structure. For more details on this please refer to the [source code](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/Structures.jl).

The following, for example, initializes a mutable structure using only the default values 

```julia
M = Geometry()
```

If one wants to modify certain parameters within the structure, one needs to call the function like: 

```julia
M = Geometry(
    ymax    =   1.0,        # [ m ]
    ymin    =   0.0,        # [ m ]
)
```

This initializes a box with a minimum and maximum depth of 0.0 [m] and 1.0 [m], respectively. The remaining parameters are defined by their default values. 

The following mutable structures including their default values are currently available: 

### Geometry
```julia 
M = Geometry(
    xmin    = 0.0,              # Minimum x-coordinate [ m ]
    xmax    = 1.0,              # Maximum x-coordinate [ m ]
    ymin    = -1.0,             # Minimum y-coordinate [ m ]
    ymax    = 0.0,              # Maximum y-coordinate [ m ]
)
```

### Physics
```julia 
P = Physics(
    g       = 9.81,             # Gravitational acceleration [ m/s² ]
    ρ₀      = 3300.0,           # Reference density [ kg/m³ ]
    k       = 4.125,            # Thermal conductivity [ W/m/K ]
    cp      = 1250.0,           # Specific heat capacity [ J/kg/K ]
    α       = 2.0e-5,           # Thermal expansion coefficient [ 1/K ]
    Q₀      = 0.0,              # Heat production rate [ W/m³ ]
    η₀      = 3.947725485e23,   # Reference viscosity [ Pa s ]
    κ       = k/ρ₀/cp,          # Thermal Diffusivity [ m²/s ]    
    ΔT      = 2500,             # Temperature difference [ K ]
    Ttop    = 273.15,           # Temperature at the top [ K ]
    Tbot    = Ttop + ΔT,        # Temperature at the bottom [ K ] 
    Ra      = 1e5,              # Rayleigh number
)
```

### Grid Spacing

```julia 
Δ = GridSpacing(
    x   =   0.0,                # Horizontal grid resolution 
    y   =   0.0,                # Vertical grid resolution
)
```

### Data Fields

```julia 
D = DataFields(
    Q       = zeros(1,1),
    T       = zeros(1,1),
    T0      = zeros(1,1),
    T_ex    = zeros(1,1),
    T_exo   = zeros(1,1),
    ρ       = zeros(1,1),
    cp      = zeros(1,1),
    vx      = zeros(1,1),
    vy      = zeros(1,1),
    Pt      = zeros(1,1),
    vxc     = zeros(1,1),
    vyc     = zeros(1,1),
    vc      = zeros(1,1),
    wt      = zeros(1,1),
    wtv     = zeros(1,1),
    ΔTtop   = zeros(1),
    ΔTbot   = zeros(1),
    Tmax    = 0.0,
    Tmin    = 0.0,
    Tmean   = 0.0,
end
)
```

If the default value for a field is used, an empty array is initialized to save memory. 

### Time Parameter

```julia
T = TimeParameter(
    const year  =   365.25*3600*24,      #   Seconds per year
    tmax        =   1000.0,              #   [ Ma ]
    Δfacc       =   0.9,                 #   Courant time factor
    Δfacd       =   0.9,                 #   Diffusion time factor
    Δ           =   0.0,                 #   Absolute time step
    Δc          =   0.0,                 #   Courant time step
    Δd          =   0.0,                 #   Diffusion time stability criterion
    itmax       =   8000,                #   Maximum iterations; 30000
)
```

Additional mutable structure will be added. For more details on the mutable structures, please refer to the [source code](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/src/Structures.jl).
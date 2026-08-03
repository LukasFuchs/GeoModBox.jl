
"""
    Geometry()

Structure containing the spatial limits of the two-dimensional model
domain.

The geometry is defined by the minimum and maximum coordinates in the
horizontal and vertical directions. The fields can be specified during
construction or modified afterward because `Geometry` is a mutable
structure.

# Fields

- `xmin`: Minimum horizontal coordinate.
- `xmax`: Maximum horizontal coordinate.
- `ymin`: Minimum vertical coordinate.
- `ymax`: Maximum vertical coordinate.

# Default values

```julia
xmin = 0.0
xmax = 1.0
ymin = -1.0
ymax = 0.0
```

# Example

```julia
M = Geometry()
```

which returns

```julia
Geometry(0.0, 1.0, -1.0, 0.0)
```

A custom model geometry can be initialized using keyword arguments:

```julia
M = Geometry(
    xmin = -50.0e3,
    xmax =  50.0e3,
    ymin =   0.0,
    ymax = 100.0e3,
)
```
"""
@kwdef mutable struct Geometry
    # Structural parameters ---
    xmin        ::Float64 = 0.0
    xmax        ::Float64 = 1.0
    ymin        ::Float64 = -1.0
    ymax        ::Float64 = 0.0
    # L           ::Float64 = 1.0
    # H           ::Float64 = 1.0
end

@kwdef mutable struct GeometryGrid
    xmin        ::Float64 = 0.0
    xmax        ::Float64 = 1.0
    ymin        ::Float64 = -1.0
    ymax        ::Float64 = 0.0
    L           ::Float64 = xmax - xmin
    H           ::Float64 = ymax - ymin
    ncx         ::Int64 = 20
    ncy         ::Int64 = 20
    nvx         ::Int64 = ncx + 1
    nvy         ::Int64 = nvy + 1 
    Δx          ::Float64 = L / ncx
    Δy          ::Float64 = H / ncy
end

"""
    Physics()

Structure containing the physical parameters required for thermal and
thermo-mechanical simulations.

The structure stores the material properties and physical constants used
throughout `GeoModBox.jl`. All quantities are specified in SI units and
can be modified during construction using keyword arguments.

# Fields

    - `g`   : Gravitational acceleration [m/s²].
    - `ρ₀`  : Reference density [kg/m³].
    - `k`   : Thermal conductivity [W/m/K].
    - `cp`  : Specific heat capacity [J/kg/K].
    - `α`   : Thermal expansion coefficient [1/K].
    - `Q₀`  : Volumetric heat-production rate [W/m³].
    - `η₀`  : Reference dynamic viscosity [Pa·s].
    - `κ`   : Thermal diffusivity [m²/s].
    - `ΔT`  : Reference temperature difference [K].
    - `Ttop`: Upper boundary temperature [K].
    - `Tbot`: Lower boundary temperature [K].
    - `Ra`  : Rayleigh number.
    - `RG`  : Universal gas constant [J/mol/K].

# Default values

```julia
g    = 9.81
ρ₀   = 3300.0
k    = 4.125
cp   = 1250.0
α    = 2.0e-5
Q₀   = 0.0
η₀   = 3.947725485e23
κ    = k / (ρ₀ * cp)
ΔT   = 2500.0
Ttop = 273.15
Tbot = Ttop + ΔT
Ra   = -9999.0
RG   = 8.314
```

# Notes

If `Ra` is negative (default), the basal Rayleigh number is computed from
the remaining physical parameters during model initialization. A positive
value may be supplied when prescribing a Rayleigh number directly.

# Example

```julia
P = Physics()

P = Physics(
    η₀ = 1.0e21,
    ΔT = 1300.0,
    Q₀ = 1.0e-6,
)
```
"""
@kwdef mutable struct Physics
    # Physical parameters --- 
    # 
    g           ::Float64 = 9.81                # Gravitational acceleration [ m/s² ]
    ρ₀          ::Float64 = 3300.0              # Reference density [ kg/m³ ]
    k           ::Float64 = 4.125               # Thermal conductivity [ W/m/K ]
    cp          ::Float64 = 1250.0              # Specific heat capacity [ J/kg/K ]
    α           ::Float64 = 2.0e-5              # Thermal expansion coefficient [ 1/K ]
    Q₀          ::Float64 = 0.0                 # Heat production rate [ W/m³ ]
    η₀          ::Float64 = 3.947725485e23      # Reference viscosity [ Pa s ]
    κ           ::Float64 = k/ρ₀/cp             # Thermal Diffusivity [ m²/s ]    
    ΔT          ::Float64 = 2500                # Temperature difference [ K ]
    Ttop        ::Float64 = 273.15              # Temperature at the top [ K ]
    Tbot        ::Float64 = Ttop + ΔT           # Temperature at the bottom [ K ] 
    Ra          ::Float64 = -9999               # Rayleigh number
    RG          ::Float64 = 8.314               # Gas Constant [ J/mol/kg ]
end

"""
    GridSpacing()

Structure containing the grid spacing of a two-dimensional Cartesian
finite-difference grid.

The structure stores the horizontal and vertical grid spacing used
throughout the numerical discretization. The values are typically
computed from the model geometry and the number of grid cells in each
coordinate direction.

# Fields

    - `x`: Horizontal grid spacing, Δx.
    - `y`: Vertical grid spacing, Δy.

# Default values

```julia
x = 0.0
y = 0.0
```

# Notes

For a regular Cartesian grid, the grid spacing is commonly computed as

```math
\\Delta x = \\frac{x_{\\max}-x_{\\min}}{N_x},
```

```math
\\Delta y = \\frac{y_{\\max}-y_{\\min}}{N_y},
```

where `Nₓ` and `Nᵧ` denote the number of cells in the horizontal and
vertical directions, respectively.

# Example

```julia
Δ = GridSpacing()

NC = (x = 100, y = 50)

Δ = GridSpacing(
    x = (M.xmax - M.xmin) / NC.x,
    y = (M.ymax - M.ymin) / NC.y,
)
```
"""
@kwdef mutable struct GridSpacing
    x           ::Float64   =   0.0
    y           ::Float64   =   0.0
end

"""
    DataFields()

Structure containing the primary model fields and auxiliary arrays used by
the numerical solvers.

`DataFields` provides a common container for thermal, mechanical, material,
and interpolation-related quantities. By default, all array fields are
initialized as `1 × 1` arrays and all scalar diagnostics are set to zero.
For an actual model, the required fields should be allocated with dimensions
consistent with the chosen finite-difference grid before the corresponding
solver is called.

# Fields

## Thermal fields

    - `Q`       : Volumetric heat-production field.
    - `Hs`      : Shear-heating field.
    - `T`       : Temperature at the cell centers.
    - `T0`      : Temperature at the previous or reference time level.
    - `T_ex`    : Extended temperature field including ghost cells.
    - `T_ex0`   : Extended temperature field at the previous or reference time level.
    - `Told_ex` : Previously stored extended temperature field.

## Viscosity fields

    - `ηc`  : Viscosity at the cell centers.
    - `η_ex`: Extended cell-centered viscosity field including ghost cells.
    - `ηv`  : Viscosity at the grid vertices.

## Density and phase fields

    - `ρ`   : Density at the cell centers.
    - `ρ_ex`: Extended cell-centered density field.
    - `p`   : Phase field at the cell centers.
    - `p_ex`: Extended cell-centered phase field.
    - `pv`  : Phase field at the grid vertices.

## Heat-capacity field

    - `cp`  : Specific heat capacity at the cell centers.

## Velocity and pressure fields

    - `vx`  : Horizontal velocity on the staggered grid.
    - `vy`  : Vertical velocity on the staggered grid.
    - `Pt`  : Pressure at the cell centers.
    - `vxc` : Horizontal velocity interpolated to the cell centers.
    - `vyc` : Vertical velocity interpolated to the cell centers.
    - `vxco`: Previously stored cell-centered horizontal velocity.
    - `vyco`: Previously stored cell-centered vertical velocity.
    - `vc`  : Cell-centered velocity magnitude.

## Interpolation weights

    - `wt`  : Generic cell-centered interpolation weights.
    - `wte` : Interpolation weights on the extended cell-centered grid.
    - `wtv` : Interpolation weights at the grid vertices.

## Thermal diagnostics

    - `ΔTtop`   : Temperature difference or temperature gradient evaluated along
                    the upper boundary.
    - `ΔTbot`   : Temperature difference or temperature gradient evaluated along
                    the lower boundary.
    - `Tmax`    : Maximum temperature.
    - `Tmin`    : Minimum temperature.
    - `Tmean`   : Mean temperature.

# Default values

All matrix fields are initialized as

```julia
zeros(1, 1)
```

all vector fields as

```julia
zeros(1)
```

and all scalar diagnostics as

```julia
0.0
```

# Notes

`DataFields` is designed as a flexible container. Individual examples and
solvers may use only a subset of the available fields. The fields required
by a specific model should therefore be allocated explicitly with dimensions
matching the corresponding staggered-grid locations.

# Example

```julia
NC = (x = 100, y = 50)
NV = (x = NC.x + 1, y = NC.y + 1)

D = DataFields(
    T    = zeros(Float64, NC.x, NC.y),
    T_ex = zeros(Float64, NC.x + 2, NC.y + 2),
    vx   = zeros(Float64, NV.x, NV.y + 1),
    vy   = zeros(Float64, NV.x + 1, NV.y),
    Pt   = zeros(Float64, NC.x, NC.y),
    ηc   = ones(Float64, NC.x, NC.y),
    ηv   = ones(Float64, NV.x, NV.y),
)
```
"""
@kwdef mutable struct DataFields
    # Data fields ---
    Q           ::Matrix{Float64}   = zeros(1,1)
    Hs          ::Matrix{Float64}   = zeros(1,1)
    T           ::Matrix{Float64}   = zeros(1,1)
    T0          ::Matrix{Float64}   = zeros(1,1)
    T_ex        ::Matrix{Float64}   = zeros(1,1)
    T_ex0       ::Matrix{Float64}   = zeros(1,1)
    Told_ex     ::Matrix{Float64}   = zeros(1,1)
    ηc          ::Matrix{Float64}   = zeros(1,1)
    η_ex        ::Matrix{Float64}   = zeros(1,1)
    ηv          ::Matrix{Float64}   = zeros(1,1)
    ρ           ::Matrix{Float64}   = zeros(1,1)
    ρ_ex        ::Matrix{Float64}   = zeros(1,1)
    p           ::Matrix{Float64}   = zeros(1,1)
    p_ex        ::Matrix{Float64}   = zeros(1,1)
    pv          ::Matrix{Float64}   = zeros(1,1)
    cp          ::Matrix{Float64}   = zeros(1,1)
    vx          ::Matrix{Float64}   = zeros(1,1)
    vy          ::Matrix{Float64}   = zeros(1,1)
    Pt          ::Matrix{Float64}   = zeros(1,1)
    vxc         ::Matrix{Float64}   = zeros(1,1)
    vyc         ::Matrix{Float64}   = zeros(1,1)
    vxco        ::Matrix{Float64}   = zeros(1,1)
    vyco        ::Matrix{Float64}   = zeros(1,1)
    vc          ::Matrix{Float64}   = zeros(1,1)
    wt          ::Matrix{Float64}   = zeros(1,1)
    wte         ::Matrix{Float64}   = zeros(1,1)
    wtv         ::Matrix{Float64}   = zeros(1,1)
    ΔTtop       ::Vector{Float64}   = zeros(1)
    ΔTbot       ::Vector{Float64}   = zeros(1)
    Tmax        ::Float64           = 0.0
    Tmin        ::Float64           = 0.0
    Tmean       ::Float64           = 0.0
end

"""
    TimeParameter()

Structure containing the parameters controlling the temporal discretization
of a simulation.

The structure stores the simulation end time, the stability factors used
to compute the Courant and diffusion time steps, the current time-step
size, and the maximum number of time steps. These parameters are shared by
the thermal, advection, and thermo-mechanical solvers.

# Fields

    - `year`    : Number of seconds in one year.
    - `tmax`    : Maximum simulation time.
    - `Δfacc`   : Safety factor for the Courant (advection) stability criterion.
    - `Δfacd`   : Safety factor for the diffusion stability criterion.
    - `Δ`       : Current time step.
    - `Δc`      : Time step determined from the Courant stability criterion.
    - `Δd`      : Time step determined from the diffusion stability criterion.
    - `itmax`   : Maximum number of iterations (time steps).

# Default values

```julia
year  = 365.25 * 24 * 3600
tmax  = 1000.0      # [Ma]
Δfacc = 0.9
Δfacd = 0.9
Δ     = 0.0
Δc    = 0.0
Δd    = 0.0
itmax = 8000
```

# Notes

The stable time step is typically chosen as

```math
\\Delta t = \\min(\\Delta t_c,\\Delta t_d),
```

where `Δt_c` is obtained from the Courant criterion and `Δt_d` from the
diffusion stability criterion. The factors `Δfacc` and `Δfacd` provide an
additional safety margin.

Although `tmax` is initialized in millions of years (Ma), it is commonly
converted to seconds before the simulation begins and subsequently
non-dimensionalized if required.

# Example

```julia
T = TimeParameter(
    tmax = 50.0,
    Δfacc = 0.8,
    Δfacd = 0.9,
    itmax = 10000,
)

T.tmax *= 1e6 * T.year
```
"""
@kwdef mutable struct TimeParameter
    # Time Parameters ---
    const year      ::Float64   =   365.25*3600*24      #   Seconds per year
    tmax            ::Float64   =   1000.0              #   [ Ma ]
    Δfacc           ::Float64   =   0.9                 #   Courant time factor
    Δfacd           ::Float64   =   0.9                 #   Diffusion time factor
    Δ               ::Float64   =   0.0                 #   Absolute time step
    Δc              ::Float64   =   0.0                 #   Courant time step
    Δd              ::Float64   =   0.0                 #   Diffusion time stability criterion
    itmax           ::Int64     =   8000                #   Maximum iterations; 30000
end
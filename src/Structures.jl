
# M = Geometry(xmin=0.0,xmax=1.0,ymin=-1.0,ymax=0.0)
# ==================================================================== #
"""
    Geometry()

Structure to initialize the geometry of the model domain.

The default values are: 

    xmin    = 0.0
    xmax    = 1.0
    ymin    = -1.0
    ymax    = 0.0

Examples
========

```julia 
julia> M = Geometry()
Geometry(0.0, 1.0, -1.0, 0.0)
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

# @kwdef mutable struct Coordinates
#     c           ::
# end

"""
    Physics()

Structure to initialize some physical constants in SI units. 

The default values are: 

    g       = 9.81                # Gravitational acceleration [ m/s² ]
    ρ₀      = 3300.0              # Reference density [ kg/m³ ]
    k       = 4.125               # Thermal conductivity [ W/m/K ]
    cp      = 1250.0              # Specific heat capacity [ J/kg/K ]
    α       = 2.0e-5              # Thermal expansion coefficient [ 1/K ]
    Q₀      = 0.0                 # Heat production rate [ W/m³ ]
    η₀      = 3.947725485e23      # Reference viscosity [ Pa s ]
    κ       = k/ρ₀/cp             # Thermal Diffusivity [ m²/s ]    
    ΔT      = 2500                # Temperature difference [ K ]
    Ttop    = 273.15              # Temperature at the top [ K ]
    Tbot    = Ttop + ΔT           # Temperature at the bottom [ K ] 
    Ra      = -9999               # Rayleigh number

The Rayleigh number is set to a negative value, which results into the calculation 
of the basal Rayleigh number based on the given physical constants. 

Examples
========

```julia
julia> P = Physics()
Physics(9.81, 3300.0, 4.125, 1250.0, 2.0e-5, 0.0, 3.947725485e23, 1.0e-6, 2500.0, 273.15, 2773.15, -9999.0)
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
end

"""
    GridSpacing()

Structure to initialize the horizontal and vertical grid spacing. 

The default values are: 

    x   =   0.0
    y   =   0.0

Examples
========

```julia
julia> Δ = GridSpacing()
GridSpacing(0.0, 0.0)
```
"""
@kwdef mutable struct GridSpacing
    x           ::Float64   =   0.0
    y           ::Float64   =   0.0
end

"""
    DataFields()

Structure to initialize some data arrays. 

The default values are: 

    Q       = zeros(1,1)
    T       = zeros(1,1)
    T0      = zeros(1,1)
    T_ex    = zeros(1,1)
    T_exo   = zeros(1,1)
    ηc      = zeros(1,1)
    η_ex    = zeros(1,1)
    ηv      = zeros(1,1)
    ρ       = zeros(1,1)
    ρ_ex    = zeros(1,1)
    cp      = zeros(1,1)
    vx      = zeros(1,1)
    vy      = zeros(1,1)
    Pt      = zeros(1,1)
    vxc     = zeros(1,1)
    vyc     = zeros(1,1)
    vc      = zeros(1,1)
    wt      = zeros(1,1)
    wte     = zeros(1,1)
    wtv     = zeros(1,1)
    ΔTtop   = zeros(1)
    ΔTbot   = zeros(1)
    Tmax    = 0.0
    Tmin    = 0.0
    Tmean   = 0.0

Examples
========
```julia
julia> D = DataFields()
DataFields([0.0;;], [0.0;;], [0.0;;], [0.0;;], [0.0;;], [0.0;;], [0.0;;], [0.0;;], [0.0;;], [0.0;;], [0.0;;], [0.0;;], [0.0;;], [0.0;;], [0.0;;], [0.0;;], [0.0;;], [0.0;;], [0.0;;], [0.0;;], [0.0], [0.0], 0.0, 0.0, 0.0)
```
"""
@kwdef mutable struct DataFields
    # Data fields ---
    Q           ::Matrix{Float64}   = zeros(1,1)
    T           ::Matrix{Float64}   = zeros(1,1)
    T0          ::Matrix{Float64}   = zeros(1,1)
    T_ex        ::Matrix{Float64}   = zeros(1,1)
    T_exo       ::Matrix{Float64}   = zeros(1,1)
    ηc          ::Matrix{Float64}   = zeros(1,1)
    η_ex        ::Matrix{Float64}   = zeros(1,1)
    ηv          ::Matrix{Float64}   = zeros(1,1)
    ρ           ::Matrix{Float64}   = zeros(1,1)
    ρ_ex        ::Matrix{Float64}   = zeros(1,1)
    cp          ::Matrix{Float64}   = zeros(1,1)
    vx          ::Matrix{Float64}   = zeros(1,1)
    vy          ::Matrix{Float64}   = zeros(1,1)
    Pt          ::Matrix{Float64}   = zeros(1,1)
    vxc         ::Matrix{Float64}   = zeros(1,1)
    vyc         ::Matrix{Float64}   = zeros(1,1)
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

Structure to initialize the time parameters. 

The default values are: 

    const year  =   365.25*3600*24      #   Seconds per year
    tmax        =   1000.0              #   [ Ma ]
    Δfacc       =   0.9                 #   Courant time factor
    Δfacd       =   0.9                 #   Diffusion time factor
    Δ           =   0.0                 #   Absolute time step
    Δc          =   0.0                 #   Courant time step
    Δd          =   0.0                 #   Diffusion time stability criterion
    itmax       =   8000                #   Maximum iterations

Examples
========

```julia
julia> T = TimeParameter
TimeParameter(3.15576e7, 1000.0, 0.9, 0.9, 0.0, 0.0, 0.0, 8000)
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

# @kwdef mutable struct Rheology

# end

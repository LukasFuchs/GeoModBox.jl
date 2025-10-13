
# M = Geometry(xmin=0.0,xmax=1.0,ymin=-1.0,ymax=0.0)
# ==================================================================== #
"""
    Geometry()

Structure to initialize the geometry of the model domain. 
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
"""
@kwdef mutable struct GridSpacing
    x           ::Float64   =   0.0
    y           ::Float64   =   0.0
end

"""
    DataFields()
"""
@kwdef mutable struct DataFields
    # Data fields ---
    Q           ::Matrix{Float64}   = zeros(1,1)
    T           ::Matrix{Float64}   = zeros(1,1)
    T0          ::Matrix{Float64}   = zeros(1,1)
    T_ex        ::Matrix{Float64}   = zeros(1,1)
    T_exo       ::Matrix{Float64}   = zeros(1,1)
    ηc          ::Matrix{Float64}   = zeros(1,1)
    ηv          ::Matrix{Float64}   = zeros(1,1)
    ρ           ::Matrix{Float64}   = zeros(1,1)
    cp          ::Matrix{Float64}   = zeros(1,1)
    vx          ::Matrix{Float64}   = zeros(1,1)
    vy          ::Matrix{Float64}   = zeros(1,1)
    Pt          ::Matrix{Float64}   = zeros(1,1)
    vxc         ::Matrix{Float64}   = zeros(1,1)
    vyc         ::Matrix{Float64}   = zeros(1,1)
    vc          ::Matrix{Float64}   = zeros(1,1)
    wt          ::Matrix{Float64}   = zeros(1,1)
    wtv         ::Matrix{Float64}   = zeros(1,1)
    ΔTtop       ::Vector{Float64}   = zeros(1)
    ΔTbot       ::Vector{Float64}   = zeros(1)
    Tmax        ::Float64           = 0.0
    Tmin        ::Float64           = 0.0
    Tmean       ::Float64           = 0.0
end

"""
    TimeParameter()
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

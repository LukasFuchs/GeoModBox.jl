"""
    Constants()

Structure to initialize the scaling constants. 

The default values are: 

    hsc = 0.0
    tsc = 0.0
    vsc = 0.0
    τsc = 0.0
    Tsc = 0.0
    Qsc = 0.0

Examples
========
```julia
julia> S = Constants()
Constants(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
```
"""
@kwdef mutable struct Constants
    hsc ::Float64 = 0.0
    tsc ::Float64 = 0.0
    vsc ::Float64 = 0.0
    τsc ::Float64 = 0.0
    Tsc ::Float64 = 0.0
    Qsc ::Float64 = 0.0
end

"""
    ScalingConstants!(M,P)

Structure to initialize and calcuate the scaling constants in SI units. 

The scaling constants are calculated by: 

    hsc =   (ymax-ymin)                 #   Length scale [ m ]
    tsc =   (ymax-ymin)^2 / κ           #   Time scale [ s ]
    vsc =   κ / (ymax-ymin)             #   Velocity scale [ m/s ]
    τsc =   (η₀ * κ)/(ymax-ymin)        #   Stress scale [ Pa ]
    Tsc =   ΔT                          #   Temperature scale [ K ]
    Qsc =   (ΔT*κ*ρ₀*cp)/(ymax-ymin)^2  #   Heat source scale [ w/m³ ]

Examples
========
```julia
julia> M = Geometry(0,1000,-1000,0)
Geometry(0.0, 1000.0, -1000.0, 0.0)

julia> P = Physics()
Physics(9.81, 3300.0, 4.125, 1250.0, 2.0e-5, 0.0, 3.947725485e23, 1.0e-6, 2500.0, 273.15, 2773.15, -9999.0)

julia> ScalingConstants!(M,P)
Constants(1000.0, 1.0e12, 9.999999999999999e-10, 3.947725485e14, 2500.0, 0.0103125)
```
"""
function ScalingConstants!(M,P)

    S       =   Constants()

    S.hsc   =   (M.ymax-M.ymin)                         #   Length scale [ m ]
    S.tsc   =   (M.ymax-M.ymin)^2 / P.κ                 #   Time scale [ s ]
    S.vsc   =   P.κ / (M.ymax-M.ymin)                   #   Velocity scale [ m/s ]
    S.τsc   =   (P.η₀ * P.κ)/(M.ymax-M.ymin)            #   Stress scale [ Pa ]
    S.Tsc   =   P.ΔT                                    #   Temperature scale [ K ]
    S.Qsc   =   (P.ΔT*P.κ*P.ρ₀*P.cp)/(M.ymax-M.ymin)^2  #   Heat source scale [ w/m³ ]

    return S
    
end

"""
    ScaleParameters!(S,M,Δ,T,P,D)

Structure to scale the parmeters using the scaling constants. 

The scaling laws are: 

    # Geometry ---
    xmin    /=  S.hsc                       # Minimum x-coordinate
    xmax    /=  S.hsc                       # Maximum x-coordinate
    ymin    /=  S.hsc                       # Minimum y-coordinate
    ymax    /=  S.hsc                       # Minimum y-coordinate
    x       /=  S.hsc                       # x-coordinate
    y       /=  S.hsc                       # y-coordinate
    # Time --- 
    tmax    /=  S.tsc                       # Maximum time
    Δc      /=  S.tsc                       # Courant time step
    Δd      /=  S.tsc                       # Diffusion time step
    Δ       /=  S.tsc                       # Total time step
    # Temperature etc. ---
    Ttop    =   (P.Ttop - 273.15)/S.Tsc     # Temperature at the top
    Tbot    =   (P.Tbot - 273.15)/S.Tsc     # Temperature at the bottom
    @. D.Q  /=  S.Qsc
    # Viscosity ---
    @. D.η  /=  P.η₀

Examples
======== 
```julia
julia> M = Geometry(0,1000,-1000,0)
Geometry(0.0, 1000.0, -1000.0, 0.0)

julia> P = Physics()
Physics(9.81, 3300.0, 4.125, 1250.0, 2.0e-5, 0.0, 3.947725485e23, 1.0e-6, 2500.0, 273.15, 2773.15, -9999.0)

julia> S = ScalingConstants!(M,P)
Constants(1000.0, 1.0e12, 9.999999999999999e-10, 3.947725485e14, 2500.0, 0.0103125)

julia> NC = (x=100,y=100,)
(x = 100, y = 100)

julia> NV      =   (
           x   =   NC.x + 1,
           y   =   NC.y + 1,
       )
(x = 101, y = 101)

julia> Δ       =   GridSpacing(
           x   =   (M.xmax - M.xmin)/NC.x,
           y   =   (M.ymax - M.ymin)/NC.y,
       )
GridSpacing(10.0, 10.0)

julia> D       =   DataFields(
           Q       =   zeros(Float64,(NC...)),
           T       =   zeros(Float64,(NC...)),
           T_ex    =   zeros(Float64,(NC.x+2,NC.y+2)),
           ρ       =   ones(Float64,(NC...)),
           cp      =   ones(Float64,(NC...)),
           vx      =   zeros(Float64,(NV.x,NV.y+1)),
           vy      =   zeros(Float64,(NV.x+1,NV.y)),    
           Pt      =   zeros(Float64,(NC...)),
           vxc     =   zeros(Float64,(NC...)),
           vyc     =   zeros(Float64,(NC...)),
           vc      =   zeros(Float64,(NC...)),
           ΔTtop   =   zeros(Float64,NC.x),
           ΔTbot   =   zeros(Float64,NC.x),
           Tmax    =   0.0,
           Tmin    =   0.0,
           Tmean   =   0.0,
       )
DataFields([0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0], [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0], [0.0;;], [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0], [0.0;;], [0.0;;], [0.0;;], [0.0;;], [1.0 1.0 … 1.0 1.0; 1.0 1.0 … 1.0 1.0; … ; 1.0 1.0 … 1.0 1.0; 1.0 1.0 … 1.0 1.0], [0.0;;], [1.0 1.0 … 1.0 1.0; 1.0 1.0 … 1.0 1.0; … ; 1.0 1.0 … 1.0 1.0; 1.0 1.0 … 1.0 1.0], [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0], [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0], [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0], [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0], [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0], [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0], [0.0;;], [0.0;;], [0.0;;], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 0.0, 0.0, 0.0)

julia> T   =   TimeParameter(
           tmax    =   1000000.0,          #   [ Ma ]
           Δfacc   =   1.0,                #   Courant time factor
           Δfacd   =   1.0,                #   Diffusion time factor
           itmax   =   8000,               #   Maximum iterations
       )
TimeParameter(3.15576e7, 1.0e6, 1.0, 1.0, 0.0, 0.0, 0.0, 8000)

julia> T.tmax      =   T.tmax*1e6*T.year    #   [ s ]
3.15576e19

julia> T.Δc        =   T.Δfacc * minimum((Δ.x,Δ.y)) / 
                           (sqrt(maximum(abs.(D.vx))^2 + maximum(abs.(D.vy))^2))
Inf

julia> T.Δd        =   T.Δfacd * (1.0 / (2.0 * P.κ *(1.0/Δ.x^2 + 1/Δ.y^2)))
2.5e7

julia> T.Δ         =   minimum([T.Δd,T.Δc])
2.5e7

ScaleParameters!(S,M,Δ,T,P,D)
```
"""
function ScaleParameters!(S,M,Δ,T,P,D)
    # Geometry ---
    M.xmin   /=  S.hsc
    M.xmax   /=  S.hsc
    M.ymin   /=  S.hsc
    M.ymax   /=  S.hsc
    Δ.x      /=  S.hsc
    Δ.y      /=  S.hsc         
    # Time --- 
    T.tmax   /=  S.tsc
    T.Δc     /=  S.tsc
    T.Δd     /=  S.tsc
    T.Δ      /=  S.tsc
    # Temperature etc. ---
    P.Ttop   =   (P.Ttop - 273.15)/S.Tsc
    P.Tbot   =   (P.Tbot - 273.15)/S.Tsc
    @. D.Q      /=  S.Qsc
    # # Viscosity ---
    # @. D.η      /=  P.η₀
end
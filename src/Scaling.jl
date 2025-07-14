"""
    Constants
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
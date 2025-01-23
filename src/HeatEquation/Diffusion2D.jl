using GeoModBox.HeatEquation.TwoD

function Advection(type,D,P,Δ,T,NC,BC)

    if type==:explicit
        ForwardEuler2Dc!(D, P.κ, Δ.x, Δ.y, T.Δ[1], P.ρ, P.cp, NC, BC)
    #elseif type==:implicit

    end

end
"""
    Itp1D_Centers2Markers!( Tm, xm, Tc, xc, Δx, xmin )
"""
function Itp1D_Centers2Markers!( Tm, xm, Tc, xc, Δx, xmin )
    # Biliniear
    for i = 1:length(xm)
        iW    = trunc(Int, (xm[i] - xmin) / Δx) + 1
        wW    = 1 - (xm[i] - xc[iW]) / Δx
        Tm[i] = wW * Tc[iW] + (1 - wW) * Tc[iW + 1]
    end
end

"""
    Itp1D_Markers2Centers!( Tc, xc, Tm, xm, dx, xmin )
"""
function Itp1D_Markers2Centers!( Tc, xc, Tm, xm, dx, xmin )
    Tc .= 0.0 
    Wc = zero(Tc)
    for i = 1:length(xm)
        iC      = trunc(Int, (xm[i] - xmin) / dx) + 1
        w       = 1 - (xm[i] - xc[iC]) / dx
        Tc[iC] += w*Tm[i]
        Wc[iC] += w
    end
    Tc ./= Wc
end
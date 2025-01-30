@doc raw"""
    VxFromVxNodes(Vx, k, p, x, y, Δ, NC, new)
"""
function VxFromVxNodes(Vx, k, p, x, y, Δ, NC, new)
    # Interpolate vx
    # From https://github.com/tduretz/M2Dpt_Julia/blob/master/Markers2D/Main_Taras_v6_Hackathon.jl
    # ---------------------------------------------------------------- #
    # Get Index of the Node ---
    i   =   Int64(round(trunc( (p.x[k] -  x.v[1])/Δ.x ) + 1))
    j   =   Int64(round(trunc( (p.y[k] - y.ce[1])/Δ.y ) + 1))
    if i<1
        i = 1
    elseif i>NC.x
        i = NC.x
    end
    if j<1
        j = 1
    elseif j> NC.y+1
        j = NC.y+1
    end
    # ---------------------------------------------------------------- #
    # Compute distances ---------------------------------------------- #
    Δxmj    =   p.x[k] -  x.v[i]
    Δymi    =   p.y[k] - y.ce[j]
    # ---------------------------------------------------------------- #
    # Compute vx velocity for the top and bottom of the cell --------- #
    vxm13   =   Vx[i,j  ] * (1-Δxmj/Δ.x) + Vx[i+1,j  ]*Δxmj/Δ.x
    vxm24   =   Vx[i,j+1] * (1-Δxmj/Δ.x) + Vx[i+1,j+1]*Δxmj/Δ.x
    # ---------------------------------------------------------------- #
    if new==1 
        if dxmj/dx>=0.5
            if i<ncx
                vxm13 += 0.5*((dxmj/dx-0.5)^2) * (Vx[i,j  ] - 2.0*Vx[i+1,j  ] + Vx[i+2,j  ])
                vxm24 += 0.5*((dxmj/dx-0.5)^2) * (Vx[i,j+1] - 2.0*Vx[i+1,j+1] + Vx[i+2,j+1])
            end
        else
            if i>1
                vxm13 += 0.5*((dxmj/dx-0.5)^2) * (Vx[i-1,j  ] - 2.0*Vx[i,j  ] + Vx[i+1,j  ])
                vxm24 += 0.5*((dxmj/dx-0.5)^2) * (Vx[i-1,j+1] - 2.0*Vx[i,j+1] + Vx[i+1,j+1])
            end
        end
    end
    # ----------------------------------------------------------------- #
    # Compute vx ------------------------------------------------------ #
    vxm     =   (1-Δymi/Δ.y) * vxm13 + (Δymi/Δ.y) * vxm24
    # ----------------------------------------------------------------- #
    return vxm
end

@doc raw"""
    VyFromVyNodes(Vy, k, p, x, y, Δ, NC, new)
"""
function VyFromVyNodes(Vy, k, p, x, y, Δ, NC, new)
    # Interpolate vy
    # From https://github.com/tduretz/M2Dpt_Julia/blob/master/Markers2D/Main_Taras_v6_Hackathon.jl
    # ---------------------------------------------------------------- #
    # Get Index of the Node ---
    i   =   Int64(round(trunc( (p.x[k] - x.ce[1])/Δ.x ) + 1))
    j   =   Int64(round(trunc( (p.y[k] -  y.v[1])/Δ.y ) + 1))
    if i<1
        i = 1
    elseif i>NC.x+1
        i = NC.x+1
    end
    if j<1
        j = 1
    elseif j>NC.y
        j = NC.y
    end
    # ---------------------------------------------------------------- #
    # Compute distances ---------------------------------------------- #
    dxmj = p.x[k] - xce[i]
    dymi = p.y[k] -  yv[j]
    # ---------------------------------------------------------------- #
    # Compute vy velocity for the left and right of the cell --------- #
    vym12 = Vy[i,j  ]*(1-dymi/dy) + Vy[i  ,j+1]*dymi/dy
    vym34 = Vy[i+1,j]*(1-dymi/dy) + Vy[i+1,j+1]*dymi/dy
    # ---------------------------------------------------------------- #
    if new==1 
        if dymi/dy>=0.5
            if j<ncy
                vym12 += 0.5*((dymi/dy-0.5)^2) * ( Vy[i,j  ] - 2.0*Vy[i,j+1  ] + Vy[i,j+2  ]);
                vym34 += 0.5*((dymi/dy-0.5)^2) * ( Vy[i+1,j] - 2.0*Vy[i+1,j+1] + Vy[i+1,j+2]);
            end      
        else
            if j>1
                vym12 += 0.5*((dymi/dy-0.5)^2) * ( Vy[i,j-1  ] - 2.0*Vy[i,j  ] + Vy[i,j+1  ]);
                vym34 += 0.5*((dymi/dy-0.5)^2) * ( Vy[i+1,j-1] - 2.0*Vy[i+1,j] + Vy[i+1,j+1]);
            end
        end
    end
    # ---------------------------------------------------------------- #
    # Compute vy ----------------------------------------------------- #
    vym = (1-dxmj/dx)*vym12 + (dxmj/dx)*vym34
    # ---------------------------------------------------------------- #
    return vym
end

@doc raw"""
    VxVyFromPrNodes(Vxp ,Vyp, k, p, xce, yce, dx, dy, ncx, ncy)
"""
function VxVyFromPrNodes(Vxp ,Vyp, k, p, x, y, Δ, NC)
    # Interpolate vx, vy
    # ---------------------------------------------------------------- #
    i   =   Int64((trunc( (p.x[k] - x.ce[1])/Δ.x ) + 1.0))
    j   =   Int64((trunc( (p.y[k] - y.ce[1])/Δ.y ) + 1.0))
    if i<1
        i = 1
    elseif i>NC.x+1
        i = NC.x+1
    end
    if j<1
        j=1
    elseif j>NC.y+1
        j = NC.y+1
    end
    # ---------------------------------------------------------------- #
    # Compute distances ---------------------------------------------- #
    Δxmj    =   p.x[k] - x.ce[i]
    Δymi    =   p.y[k] - y.ce[j]
    # ---------------------------------------------------------------- #
    # Compute weights ------------------------------------------------ #
    wtmij   =   (1.0-Δxmj/Δ.x)*(1.0-Δymi/Δ.y)
    wtmi1j  =   (1.0-Δxmj/Δ.x)*(    Δymi/Δ.y)    
    wtmij1  =   (    Δxmj/Δ.x)*(1.0-Δymi/Δ.y)
    wtmi1j1 =   (    Δxmj/Δ.x)*(    Δymi/Δ.y)
    # ---------------------------------------------------------------- #
    # Compute vx, vy velocity ---------------------------------------- #
    vxm = Vxp[i,j]*wtmij + Vxp[i,j+1]*wtmi1j + Vxp[i+1,j]*wtmij1 + Vxp[i+1,j+1]*wtmi1j1
    vym = Vyp[i,j]*wtmij + Vyp[i,j+1]*wtmi1j + Vyp[i+1,j]*wtmij1 + Vyp[i+1,j+1]*wtmi1j1
    # ---------------------------------------------------------------- #
    return vxm, vym
end
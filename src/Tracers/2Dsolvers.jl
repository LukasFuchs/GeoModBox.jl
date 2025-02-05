mutable struct TMarkers
    x       ::  Array{Float64,1}
    y       ::  Array{Float64,1}
    T       ::  Array{Float64,1}
    mpc     ::  Array{Float64,2}
    phase   ::  Array{Int64,1}
end

@doc raw"""
    IniTracer2D(nmx,nmy,Δ,M,NC,noise)
"""
function IniTracer2D(Aparam,nmx,nmy,Δ,M,NC,noise)
    
    nmark   =   nmx*nmy*NC.x*NC.y

    # Initialise markers ---
    Δxm, Δym  =   Δ.x/nmx, Δ.y/nmy 
    
    xm  =   LinRange(M.xmin+Δxm/2.0, M.xmax-Δxm/2.0, NC.x*nmx)
    ym  =   LinRange(M.ymin+Δym/2.0, M.ymax-Δym/2.0, NC.y*nmy)        
    
    (xmi,ymi) = ([x for x=xm,y=ym], [y for x=xm,y=ym])
    
    # Over allocate markers
    xm  =   vec(xmi)
    ym  =   vec(ymi)    
    phm =   zeros(Int64,   size(xm))
    mpc =   zeros(Float64,(NC.x,NC.y))
    if Aparam==:thermal
        Tm  =   zeros(Float64, size(xm))
        Ma  =   TMarkers( xm, ym, Tm, mpc, phm )
    end

    ## define phase ---
    #for k=1:nmark
    #    if (Ma.x[k]<Ma.y[k]) 
    #        Ma.phase[k] = 1
    #    end
    #end

    # add noise ---
    if noise==1
        for k=1:nmark
            Ma.x[k] += (rand()-0.5)*Δxm
            Ma.y[k] += (rand()-0.5)*Δym
        end
    end

    return Ma
end

@doc raw"""
    CountMPC()
"""
function CountMPC(Ma,nmark,x,y,Δ)
    for k = 1:nmark
        if (Ma.phase[k]>=0)
            # Get the column ---
            dstx    =   Ma.x[k] - x.c[1]
            i       =   Int64(round(ceil( (dstx/Δ.x) + 0.5)))
            # Get the line ---
            dsty    =   Ma.y[k] - y.c[1]
            j       =   Int64(round(ceil( (dsty/Δ.y) + 0.5)))
            # Increment cell count ---
            Ma.mpc[i,j]     +=  1.0
        end
    end
    return Ma
end

@doc raw"""
    VxFromVxNodes(Vx, k, Ma, x, y, Δ, NC, new)
"""
function VxFromVxNodes(Vx, k, Ma, x, y, Δ, NC, new)
    # Interpolate vx
    # From https://github.com/tduretz/M2Dpt_Julia/blob/master/Markers2D/Main_Taras_v6_Hackathon.jl
    # ---------------------------------------------------------------- #
    # Get Index of the Node ---
    i   =   Int64(round(trunc( (Ma.x[k] -  x.v[1])/Δ.x ) + 1))
    j   =   Int64(round(trunc( (Ma.y[k] - y.ce[1])/Δ.y ) + 1))
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
    Δxmj    =   Ma.x[k] -  x.v[i]
    Δymi    =   Ma.y[k] - y.ce[j]
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
    VyFromVyNodes(Vy, k, Ma, x, y, Δ, NC, new)
"""
function VyFromVyNodes(Vy, k, Ma, x, y, Δ, NC, new)
    # Interpolate vy
    # From https://github.com/tduretz/M2Dpt_Julia/blob/master/Markers2D/Main_Taras_v6_Hackathon.jl
    # ---------------------------------------------------------------- #
    # Get Index of the Node ---
    i   =   Int64(round(trunc( (Ma.x[k] - x.ce[1])/Δ.x ) + 1))
    j   =   Int64(round(trunc( (Ma.y[k] -  y.v[1])/Δ.y ) + 1))
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
    dxmj = Ma.x[k] - xce[i]
    dymi = Ma.y[k] -  yv[j]
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
    VxVyFromPrNodes(Vxp ,Vyp, k, Ma, x, y, Δ, NC )
"""
function VxVyFromPrNodes(Vxp ,Vyp, k, Ma, x, y, Δ, NC)
    # Interpolate vx, vy
    # ---------------------------------------------------------------- #
    i   =   Int64((trunc( (Ma.x[k] - x.ce[1])/Δ.x ) + 1.0))
    j   =   Int64((trunc( (Ma.y[k] - y.ce[1])/Δ.y ) + 1.0))
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
    Δxmj    =   Ma.x[k] - x.ce[i]
    Δymi    =   Ma.y[k] - y.ce[j]
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

@doc raw"""
    FromCtoM(Prop, Ma, x, y, Δ, NC)
"""
function FromCtoM(Prop, k, Ma, x, y, Δ, NC)
     # Interpolate Property from Centroids to Marker
    # ---------------------------------------------------------------- #
    i   =   Int64((trunc( (Ma.x[k] - x.ce[1])/Δ.x ) + 1.0))
    j   =   Int64((trunc( (Ma.y[k] - y.ce[1])/Δ.y ) + 1.0))
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
    Δxmj    =   Ma.x[k] - x.ce[i]
    Δymi    =   Ma.y[k] - y.ce[j]
    # ---------------------------------------------------------------- #
    # Compute weights ------------------------------------------------ #
    wtmij   =   (1.0-Δxmj/Δ.x)*(1.0-Δymi/Δ.y)
    wtmi1j  =   (1.0-Δxmj/Δ.x)*(    Δymi/Δ.y)    
    wtmij1  =   (    Δxmj/Δ.x)*(1.0-Δymi/Δ.y)
    wtmi1j1 =   (    Δxmj/Δ.x)*(    Δymi/Δ.y)
    # ---------------------------------------------------------------- #
    # Compute Marker Property ---------------------------------------- #
    Propm   =   Prop[i,j]*wtmij + 
                    Prop[i,j+1]*wtmi1j + 
                    Prop[i+1,j]*wtmij1 + 
                    Prop[i+1,j+1]*wtmi1j1
    # ---------------------------------------------------------------- #
    return Propm
end

@doc raw"""
    Markers2Cells(Ma,nmark,PG,weight,x,y,Δ,param)
"""
function Markers2Cells(Ma,nmark,PG,weight,x,y,Δ,param)
    PG      .*=     0.0
    weight  .*=     0.0
    if param==:thermal
        PM  =       copy(Ma.T)
    end
    #chunks = Iterators.partition(1:nmark, nmark ÷ nthreads())
    #@sync for chunk in chunks
    #    @spawn begin
    #        tid = threadid()
    #        fill!(phase_th[tid], 0)
    #        fill!(weight_th[tid], 0)
            #for k in chunk
            for k = 1:nmark
                # Get the column:
                dstx = Ma.x[k] - x.c[1]
                i = ceil(Int, dstx / Δ.x + 0.5)
                # Get the line:
                dsty = Ma.y[k] - y.c[1]
                j = ceil(Int, dsty /  Δ.y + 0.5)
                # Relative distances
                Δxm = 2.0 * abs(x.c[i] - Ma.x[k])
                Δym = 2.0 * abs(y.c[j] - Ma.y[k])
                # Increment cell counts
                area = (1.0 - Δxm / Δ.x) * (1.0 - Δym / Δ.y)
                PG[i, j] += PM[k] * area
                weight[i, j] += area
            end
    #    end
    #end

    #phase  .= reduce(+, phase_th)
    #weight .= reduce(+, weight_th)
    PG ./= weight

    return
end

@doc raw"""
    FromCtoM(Prop, Ma, x, y, Δ, NC)
"""
function AdvectTracer2D()

end
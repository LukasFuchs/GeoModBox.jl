using Plots
using ExtendableSparse, GeoModBox.InitialCondition

function main()

    # Geometry ---------------------------------------------------------- #
    M       =   (
        xmin    =   0.0,
        xmax    =   500.0e3,    # [ m ]
        ymin    =   -500.0e3,   # [ m ]
        ymax    =   0.0,
    )
    # -------------------------------------------------------------------- #
    # Grid --------------------------------------------------------------- #
    NC      =   (
        x   =   10, 
        y   =   10,
    )
    NV      =   (
        x   =   NC.x + 1,
        y   =   NC.y + 1,
    )
    Δ       =   (
        x   =   (M.xmax - M.xmin)/NC.x,
        y   =   (M.ymax - M.ymin)/NC.y,
    )
    x       =   (
        c   =   LinRange(M.xmin+Δ.x/2,M.xmax-Δ.x/2,NC.x),
        ce  =   LinRange(M.xmin - Δ.x/2.0, M.xmax + Δ.x/2.0, NC.x+2),
        v   =   LinRange(M.xmin,M.xmax,NV.x),
    )
    y       =   (
        c   =   LinRange(M.ymin+Δ.y/2,M.ymax-Δ.y/2,NC.y),
        ce  =   LinRange(M.ymin - Δ.x/2.0, M.ymax + Δ.x/2.0, NC.y+2),
        v   =   LinRange(M.ymin,M.ymax,NV.y),
    )
    x1      =   (
        c2d     =   x.c .+ 0*y.c',
        v2d     =   x.v .+ 0*y.v', 
        vx2d    =   x.v .+ 0*y.ce',
        vy2d    =   x.ce .+ 0*y.v',
    )
    x   =   merge(x,x1)
    y1      =   (
        c2d     =   0*x.c .+ y.c',
        v2d     =   0*x.v .+ y.v',
        vx2d    =   0*x.v .+ y.ce',
        vy2d    =   0*x.ce .+ y.v',
    )
    y   =   merge(y,y1)
    # -------------------------------------------------------------------- #
    # Physics ------------------------------------------------------------ #
    g       =   9.81

    η₀      =   1.0e21
    #η₁      =   1.0e21
    #η       =   [η₀,η₁]

    ρ₀      =   3200.0          #   Background density
    ρ₁      =   3300.0          #   Block density
    ρ       =   [ρ₀,ρ₁]    
    
    phase   =   [0,1]
    # ------------------------------------------------------------------- #
    # Allocation -------------------------------------------------------- #
    D   =   (
        vx  =   zeros(NV.x,NC.y+2),
        vy  =   zeros(NC.x+2,NV.y),
        Pt  =   zeros(NC...),
        ρ   =   zeros(NC...),
        vxc =   zeros(NC...),
        vyc =   zeros(NC...),
        vc  =   zeros(NC...),
    )
    # ------------------------------------------------------------------- #
    # Initial Condition ------------------------------------------------- #
    type    =:block
    IniPhase!(type,D,M,x,y,NC;ρ)
    p   =   heatmap(x.c./1e3,y.c./1e3,D.ρ',color=:inferno,
                xlabel="x[km]",ylabel="y[km]",
                aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3), 
                ylims=(M.ymin/1e3, M.ymax/1e3),colorbar=true,
                layout=(3,1),subplot=1)
    display(p)
    # ------------------------------------------------------------------- #
    # Numbering, without ghost nodes! ---
    off    = [  NV.x*NC.y,                          # vx
                NV.x*NC.y + NC.x*NV.y,              # vy
                NV.x*NC.y + NC.x*NV.y + NC.x*NC.y]  # Pt

    Num    =    (
        Vx  =   reshape(1:NV.x*NC.y, NV.x, NC.y), 
        Vy  =   reshape(off[1]+1:off[1]+NC.x*NV.y, NC.x, NV.y), 
        Pt  =   reshape(off[2]+1:off[2]+NC.x*NC.y,NC...),
    )

    rhs     =   zeros(maximum(Num.Pt))
    χ       =   zeros(maximum(Num.Pt))

    VBC     =   (
        type    =   (E=:freeslip,W=:freeslip,S=:freeslip,N=:freeslip),
        val     =   (E=zeros(NV.y),W=zeros(NV.y),S=zeros(NV.x),N=zeros(NV.x)),
    )

    K       =   Assemblyc(NC, NV, Δ, η₀, VBC, Num)

    rhs     =   updaterhs( NC, NV, Δ, η₀, D.ρ, g, VBC, Num, rhs )
    
    # Solution --- 
    χ       =   K \ rhs

    @show size(χ)

    D.vx[:,2:end-1]     .=  χ[Num.Vx]
    D.vy[2:end-1,:]     .=  χ[Num.Vy]
    D.Pt                .=  χ[Num.Pt]

    # Get the velocity on the centroids ---
    for i = 1:NC.x
        for j = 1:NC.y
            D.vxc[i,j]  = (D.vx[i,j+1] + D.vx[i+1,j+1])/2
            D.vyc[i,j]  = (D.vy[i+1,j] + D.vy[i+1,j+1])/2
        end
    end
    @. D.vc        = sqrt(D.vxc^2 + D.vyc^2)

    @show(minimum(D.vc))
    @show(maximum(D.vc))

    quiver!(p,x.c2d./1e3,y.c2d./1e3,
            quiver=(D.vx,D.vyc),        
            color="white",layout=(3,1),subplot=1)
    heatmap!(p,x.c./1e3,y.c./1e3,D.vxc',
                aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),
                ylims=(M.ymin/1e3, M.ymax/1e3),
                layout=(3,1),subplot=2)
    heatmap!(p,x.c./1e3,y.c./1e3,D.vyc',
                aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),
                ylims=(M.ymin/1e3, M.ymax/1e3),
                layout=(3,1),subplot=3)
    display(p)

end

function Assemblyc(NC, NV, Δ, η, BC, Num)

    # Linear system of equation ---
    ndof    =   maximum(Num.Pt)
    K       =   ExtendableSparseMatrix(ndof,ndof)
    dx,dy   =   Δ.x, Δ.y

    # x momentum equation ----------------------------------------------- #
    for i = 1:NV.x, j = 1:NC.y
        # Equation Number ---
        ii  =   Num.Vx[i,j] 
        if i == 1 || i == NV.x
            # East and West boundary ---
            # Free Slip && No Slip: vₓ = 0 
            K[ii,ii]    =   1.0
        else
            # Stencil, vₓ ---
            iS  =   ii - NV.x
            iW  =   ii - 1 
            iC  =   ii 
            iE  =   ii + 1
            iN  =   ii + NV.x
            # Pressure ---
            iPE =   Num.Pt[i,j]
            iPW =   Num.Pt[i-1,j]
            # ---
            inS     =   j==1    ? false  : true  
            FSS     =   (j==1    && BC.type.S==:freeslip) ? 1. : 0.
            NSS     =   (j==1    && BC.type.S==:noslip) ? 1. : 0.
            inN     =   j==NC.y ? false  : true   
            FSN     =   (j==NC.y && BC.type.N==:freeslip) ? 1. : 0.
            NSN     =   (j==NC.y && BC.type.N==:noslip) ? 1. : 0.
            if inS K[ii,iS] =   η / dy^2 end
            K[ii,iW]    =   η / dx^2
            K[ii,iC]    =   - 2.0 * η / dx^2 - (2.0 - FSS - FSN + NSS + NSN) * η / dy^2
            K[ii,iE]    =   η / dx^2
            if inN K[ii,iN] =   η / dy^2 end
            K[ii,iPW]   = 1 / dx
            K[ii,iPE]   = -1 / dx
        end
    end
    # ------------------------------------------------------------------- #

    # y momentum equation ----------------------------------------------- #
    for i = 1:NC.x, j = 1:NV.y
        # Equation Number ---
        ii  =   Num.Vy[i,j] 
        if j == 1 || j == NV.y
            # East and West boundary ---
            # Free Slip && No Slip: vₓ = 0 
            K[ii,ii]    =   1.0
        else
            # Stencil, vₓ ---
            iS  =   ii - NC.x
            iW  =   ii - 1 
            iC  =   ii 
            iE  =   ii + 1
            iN  =   ii + NC.x
            # Pressure ---
            iPS =   Num.Pt[i,j-1]   
            iPN =   Num.Pt[i,j]
            # ---
            inW     =   i==1    ? false  : true  
            FSW     =   (i==1    && BC.type.W==:freeslip) ? 1. : 0.
            NSW     =   (i==1    && BC.type.W==:noslip) ? 1. : 0.
            inE     =   i==NC.x ? false  : true   
            FSE     =   (i==NC.x && BC.type.E==:freeslip) ? 1. : 0.
            NSE     =   (i==NC.x && BC.type.E==:noslip) ? 1. : 0.            
            K[ii,iS]    =   η / dy^2
            if inW K[ii,iW] =  η / dx^2 end
            K[ii,iC]    =   - (2.0 - FSW - FSE + NSW + NSE) * η / dx^2 - 2.0 * η / dy^2
            if inE K[ii,iE] = η / dx^2 end
            K[ii,iN]    =   η / dy^2
            K[ii,iPS]   =   1 / dy
            K[ii,iPN]   =   - 1 / dy
        end
    end
    # ------------------------------------------------------------------- #

    # Continuity equation ----------------------------------------------- #
    for i = 1:NC.x, j = 1:NC.y
        # Equation number ---
        ii  =   Num.Pt[i,j]
        # Stencil ---
        iW  =   Num.Vx[i,j]
        iE  =   Num.Vx[i+1,j]
        iS  =   Num.Vy[i,j]
        iN  =   Num.Vy[i,j+1]
        # Linear system coefficients
        K[ii,iW]    =   -1 / dx
        K[ii,iE]    =   1 / dx
        K[ii,iS]    =   -1 / dy 
        K[ii,iN]    =   1 / dy
    end

    return flush!(K)
end

function updaterhs( NC, NV, Δ, η, ρ, g, BC, Num, rhs)

    # x momentum equation ----------------------------------------------- #
    for i = 1:NV.x, j = 1:NC.y
        # Equation Number ---
        ii  =   Num.Vx[i,j] 
        if i == 1 || i == NV.x
            # East and West boundary ---
            # Free Slip && No Slip: vₓ = 0 
            rhs[ii]     =   0.0
        else            
            # ---
            NSS     =   (j==1    && BC.type.S==:noslip) ? 1. : 0.
            NSN     =   (j==NC.y && BC.type.N==:noslip) ? 1. : 0.            
            # ---
            rhs[ii] +=  -2.0 * η * BC.val.S[i] / Δ.y^2 * NSS -
                        2.0 * η * BC.val.N[i] / Δ.y^2 * NSN
        end
    end
    # y momentum equation ----------------------------------------------- #
    for i = 1:NC.x, j = 1:NV.y
        # Equation Number ---
        ii  =   Num.Vy[i,j] 
        if j == 1 || j == NV.y
            # North and South boundary ---
            # Free Slip && No Slip: vy = 0 
            rhs[ii]     =   0.0
        else            
            # ---
            NSW     =   (i==1    && BC.type.W==:noslip) ? 1. : 0.
            NSE     =   (i==NC.x && BC.type.E==:noslip) ? 1. : 0.            
            # ---
            rhs[ii] +=  g * ((ρ[i,j] + ρ[i,j-1]) / 2.0) - 
                        2.0 * η * BC.val.W[i] / Δ.x^2 * NSW -
                        2.0 * η * BC.val.E[i] / Δ.x^2 * NSE
        end
    end

    return rhs
end


main()
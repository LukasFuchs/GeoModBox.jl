using ExtendableSparse

function Assemblyc(η, BC, Num, NV, Nc, Δ)

    # Linear system of equation ---
    ndof    =   maximum(Num.Pt)
    K       =   ExtendableSparse(ndof,ndof)
    dx,dy   =   Δ.x, Δ.y

    # x momentum equation ----------------------------------------------- #
    for i = 1:NV.x, j = 1:NC.y
        
        # Equation Number ---
        ii  =   Num.Vx[i,j] 
        # Stencil, vₓ ---
        iS  =   ii - NV.x
        iW  =   ii - 1 
        iC  =   ii 
        iW  =   ii + 1
        iN  =   ii + NV.x
        # Pressure ---
        iPE =   Num.Pt[i,j]
        iPW =   Num.Pt[i-1,j]

    end

end
using ExtendableSparse

function Poisson2Dc!(D,NC,P,BC,Δ,K,rhs,Num)
# Function to solve 2D heat diffusion equation using the explicit finite
# difference scheme
# [Q] = W/m^3
# ----------------------------------------------------------------------- #
    
    # a       =   1.0 / Δ.x[1]^2.0
    # b       =   1.0 / Δ.y[1]^2.0

    #  --------------------------------------------- #
    rhs     .=   - reshape(D.Q, NC.xc*NC.yc, 1) ./ P.k

    for i=1:NC.xc, j=1:NC.yc
        # Equation number
        ii = Num.T[i,j]
        # Stencil
        iS = ii - NC.xc
        iW = ii - 1         
        iC = ii        
        iE = ii + 1
        iN = ii + NC.xc
        # Boundaries
        # West boundary ---
        inW    =  i==1    ? false  : true   
        DirW   = (i==1    && BC.type.W==:Dirichlet) ? 1. : 0.
        NeuW   = (i==1    && BC.type.W==:Neumann  ) ? 1. : 0.
        # East boundary ---
        inE    =  i==NC.xc ? false  : true   
        DirE   = (i==NC.xc && BC.type.E==    :Dirichlet) ? 1. : 0.
        NeuE   = (i==NC.xc && BC.type.E==:Neumann  ) ? 1. : 0.
        # South boundary ---
        inS    =  j==1    ? false  : true  
        DirS   = (j==1    && BC.type.S==:Dirichlet) ? 1. : 0.
        NeuS   = (j==1    && BC.type.S==:Neumann  ) ? 1. : 0.
        # North boundary ---
        inN    =  j==NC.yc ? false  : true   
        DirN   = (j==NC.yc && BC.type.N==:Dirichlet) ? 1. : 0.
        NeuN   = (j==NC.yc && BC.type.N==:Neumann  ) ? 1. : 0.
        # Linear system coefficients
        if inS K[ii,iS] =  1.0 / Δ.y[1]^2.0 end
        if inW K[ii,iW] =  1.0 / Δ.x[1]^2.0 end
        K[ii,iC] = - ( (2.0 + DirW + DirE - NeuW - NeuE) / Δ.x[1]^2.0 + (2.0 + DirS + DirN - NeuS - NeuN) / Δ.y[1]^2.0 )
        if inE K[ii,iE] = 1.0 / Δ.x[1]^2.0 end
        if inN K[ii,iN] = 1.0 / Δ.y[1]^2.0 end
        
        # Update right hand side 
        rhs[ii]     += - 2.0 * BC.val.W * DirW / Δ.x[1]^2.0 
                        - 2.0 * BC.val.E * DirE / Δ.x[1]^2.0 
                        - 2.0 * BC.val.S * DirS / Δ.y[1]^2.0 
                        - 2.0 * BC.val.N * DirN / Δ.y[1]^2.0
                        + BC.val.W * Δ.x[1] * NeuW / Δ.x[1]^2.0
                        - BC.val.E * Δ.x[1] * NeuE / Δ.x[1]^2.0
                        + BC.val.S * Δ.y[1] * NeuS / Δ.y[1]^2.0
                        - BC.val.N * Δ.y[1] * NeuN / Δ.y[1]^2.0
    end

    D.T[:]  .=   K \ rhs[:]
    
end

function Poisson2D!( T, Q, kx, ky, Δx, Δy, NC, BC, K, rhs, Num ) 
    #  --------------------------------------------- #
    rhs     .=   - reshape(Q, NC.xc*NC.yc, 1)

    for i=1:NC.xc, j=1:NC.yc
        # Equation number
        ii = Num.T[i,j]
        # Stencil
        iS = ii - NC.xc
        iW = ii - 1         
        iC = ii        
        iE = ii + 1
        iN = ii + NC.xc
        # Boundaries
        # West boundary ---
        inW    =  i==1    ? false  : true   
        DirW   = (i==1    && BC.type.W==:Dirichlet) ? 1. : 0.
        NeuW   = (i==1    && BC.type.W==:Neumann  ) ? 1. : 0.
        # East boundary ---
        inE    =  i==NC.xc ? false  : true   
        DirE   = (i==NC.xc && BC.type.E==    :Dirichlet) ? 1. : 0.
        NeuE   = (i==NC.xc && BC.type.E==:Neumann  ) ? 1. : 0.
        # South boundary ---
        inS    =  j==1    ? false  : true  
        DirS   = (j==1    && BC.type.S==:Dirichlet) ? 1. : 0.
        NeuS   = (j==1    && BC.type.S==:Neumann  ) ? 1. : 0.
        # North boundary ---
        inN    =  j==NC.yc ? false  : true   
        DirN   = (j==NC.yc && BC.type.N==:Dirichlet) ? 1. : 0.
        NeuN   = (j==NC.yc && BC.type.N==:Neumann  ) ? 1. : 0.
        # Linear system coefficients
        if inS K[ii,iS] = ky[i,j] / Δy^2.0 end
        if inW K[ii,iW] = kx[i,j] / Δx^2.0 end
        K[ii,iC] = - ( kx[i,j] * (1.0 + DirW - NeuW) + kx[i+1,j] * (1.0 + DirE - NeuE) ) / Δx^2.0 - ( ky[i,j] * (1.0 + DirS - NeuS) + ky[i,j+1] * (1.0 + DirN - NeuN) ) / Δy^2.0
        if inE K[ii,iE] = kx[i+1,j] / Δx^2.0 end
        if inN K[ii,iN] = ky[i,j+1] / Δy^2.0 end
        
        # Update right hand side 
        rhs[ii]     += - 2.0 * BC.val.W[j] * kx[i,j] * DirW / Δx^2.0 
                        - 2.0 * BC.val.E[j] * kx[i+1,j] * DirE / Δx^2.0
                        - 2.0 * BC.val.S[i] * ky[i,j] * DirS / Δy^2.0 
                        - 2.0 * BC.val.N[i] * ky[i,j+1] * DirN / Δy^2.0
                        + BC.val.W[j] * Δx * kx[i,j] * NeuW / Δx^2.0
                        - BC.val.E[j] * Δx * kx[i+1,j] * NeuE / Δx^2.0
                        + BC.val.S[i] * Δy * ky[i,j] * NeuS / Δy^2.0
                        - BC.val.N[i] * Δy * ky[i,j+1] * NeuN / Δy^2.0
    end

    T[:]  .=   K \ rhs[:]
end
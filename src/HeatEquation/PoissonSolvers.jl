using ExtendableSparse

function Poisson!(D,NC,P,BC,Δ,K,rhs,Num)
# Function to solve 2D heat diffusion equation using the explicit finite
# difference scheme
# [Q] = W/m^3
# ----------------------------------------------------------------------- #
    
    # a       =   1.0 / Δ.x[1]^2.0
    # b       =   1.0 / Δ.y[1]^2.0

    #  --------------------------------------------- #
    rhs     .=   - reshape(D.Q, NC.x*NC.y, 1) ./ P.k

    for i=1:NC.x, j=1:NC.y
        # Equation number
        ii = Num.T[i,j]
        # Stencil
        iS = ii - NC.x
        iW = ii - 1         
        iC = ii        
        iE = ii + 1
        iN = ii + NC.x
        # Boundaries
        # West boundary ---
        inW    =  i==1    ? false  : true   
        DirW   = (i==1    && BC.type.W==:Dirichlet) ? 1. : 0.
        NeuW   = (i==1    && BC.type.W==:Neumann  ) ? 1. : 0.
        # East boundary ---
        inE    =  i==NC.x ? false  : true   
        DirE   = (i==NC.x && BC.type.E==    :Dirichlet) ? 1. : 0.
        NeuE   = (i==NC.x && BC.type.E==:Neumann  ) ? 1. : 0.
        # South boundary ---
        inS    =  j==1    ? false  : true  
        DirS   = (j==1    && BC.type.S==:Dirichlet) ? 1. : 0.
        NeuS   = (j==1    && BC.type.S==:Neumann  ) ? 1. : 0.
        # North boundary ---
        inN    =  j==NC.y ? false  : true   
        DirN   = (j==NC.y && BC.type.N==:Dirichlet) ? 1. : 0.
        NeuN   = (j==NC.y && BC.type.N==:Neumann  ) ? 1. : 0.
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
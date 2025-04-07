using ExtendableSparse

@doc raw"""
    ForwardEuler2Dc
"""
function ForwardEuler2Dc!(D, κ, Δx, Δy, Δt, ρ, cp, NC, BC)
    # Function to solve 2D heat diffusion equation using the explicit finite
    # difference scheme
    # Q - Waermeproduktionsrate pro Volumen [W/m^3]
    # ------------------------------------------------------------------- #
    
    sx      = κ * Δt / Δx^2
    sz      = κ * Δt / Δy^2

    # Temperature at the ghost nodes ------------------------------------ #
    # West boundary ---
    D.T_ex[1,2:end-1]   .= (BC.type.W==:Dirichlet) .* (2 .* BC.val.W .- D.T_ex[2,2:end-1]) + 
                            (BC.type.W==:Neumann) .* (D.T_ex[2,2:end-1] .- BC.val.W .* Δx)
    # East boundary ---
    D.T_ex[end,2:end-1] .= (BC.type.E==:Dirichlet) .* (2 .* BC.val.E .- D.T_ex[end-1,2:end-1]) + 
                            (BC.type.E==:Neumann) .* (D.T_ex[end-1,2:end-1] .+ BC.val.E .* Δx)
    # South boundary --- 
    D.T_ex[2:end-1,1]   .= (BC.type.S==:Dirichlet) .* (2 .* BC.val.S .- D.T_ex[2:end-1,2]) + 
                            (BC.type.S==:Neumann) .* (D.T_ex[2:end-1,2] .- BC.val.S .* Δy)
    # Northern boundary ---
    D.T_ex[2:end-1,end] .= (BC.type.N==:Dirichlet) .* (2 .* BC.val.N .- D.T_ex[2:end-1,end-1]) + 
                            (BC.type.N==:Neumann) .* (D.T_ex[2:end-1,end-1] .+ BC.val.N .* Δy)
    # ------------------------------------------------------------------- #
    # Loop over internal nodes ------------------------------------------ #
    for i = 1:NC.x
        for j = 1:NC.y
            i1 = i+1
            j1 = j+1
            D.T[i,j] = D.T_ex[i1,j1] + 
                sx * (D.T_ex[i1-1,j1] - 2 * D.T_ex[i1,j1] + D.T_ex[i1+1,j1]) + 
                sz * (D.T_ex[i1,j1-1] - 2 * D.T_ex[i1,j1] + D.T_ex[i1,j1+1]) + 
                D.Q[i,j] * Δt / ρ[i,j] / cp
        end
    end
    # ------------------------------------------------------------------- #
    # Update extended temperature --------------------------------------- #
    D.T_ex[2:end-1,2:end-1]     .=  D.T
    # ------------------------------------------------------------------- #    
end

@doc raw"""
    ComputeResiduals2D
"""
function ComputeResiduals2D!(R, T, T_ex, T0, ∂T, q, ρ, Cp, k, BC, Δ, Δt)
    @. T_ex[2:end-1,2:end-1] = T 
    @. T_ex[  1,2:end-1] = (BC.type.W==:Dirichlet) * (2*BC.val.W - T_ex[    2,2:end-1]) + (BC.type.W==:Neumann) * (T_ex[    2,2:end-1] - Δ.x/k.x[  1,:]*BC.val.W)
    @. T_ex[end,2:end-1] = (BC.type.E==:Dirichlet) * (2*BC.val.E - T_ex[end-1,2:end-1]) + (BC.type.E==:Neumann) * (T_ex[end-1,2:end-1] + Δ.x/k.x[end,:]*BC.val.E)
    @. T_ex[2:end-1,  1] = (BC.type.S==:Dirichlet) * (2*BC.val.S - T_ex[2:end-1,    2]) + (BC.type.S==:Neumann) * (T_ex[2:end-1,    2] - Δ.y/k.y[:,  1]*BC.val.S)
    @. T_ex[2:end-1,end] = (BC.type.N==:Dirichlet) * (2*BC.val.N - T_ex[2:end-1,end-1]) + (BC.type.N==:Neumann) * (T_ex[2:end-1,end-1] - Δ.y/k.y[:,end]*BC.val.N)
    @. ∂T.∂x = (T_ex[2:end,2:end-1] - T_ex[1:end-1,2:end-1])/Δ.x
    @. ∂T.∂y = (T_ex[2:end-1,2:end] - T_ex[2:end-1,1:end-1])/Δ.y
    @. q.x   = -k.x * ∂T.∂x
    @. q.y   = -k.y * ∂T.∂y
    @. R     = ρ*Cp*(T - T0)/Δt + (q.x[2:end,:] - q.x[1:end-1,:])/Δ.x + (q.y[:,2:end] - q.y[:,1:end-1])/Δ.y  
end

@doc raw"""
    AssembleMatrix2D
"""
function AssembleMatrix2D(rho, cp, k, BC, Num, nc, Δ, Δt)
    # Linear system of equation
    ndof   = maximum(Num.T)
    K      = ExtendableSparseMatrix(ndof, ndof)
    dx, dy = Δ.x, Δ.y
    #############################
    #       Heat equation       #
    #############################
    for i=1:nc.x
        for j=1:nc.y
            # Equation number
            ii = Num.T[i,j]
            # Stencil
            iS = ii - nc.x
            iW = ii - 1
            iC = ii
            iE = ii + 1
            iN = ii + nc.x
            # Boundaries
            inW    =  i==1    ? false  : true   
            DirW   = (i==1    && BC.type.W==:Dirichlet) ? 1. : 0.
            NeuW   = (i==1    && BC.type.W==:Neumann  ) ? 1. : 0.
            inE    =  i==nc.x ? false  : true   
            DirE   = (i==nc.x && BC.type.E==:Dirichlet) ? 1. : 0.
            NeuE   = (i==nc.x && BC.type.E==:Neumann  ) ? 1. : 0.
            inS    =  j==1    ? false  : true  
            DirS   = (j==1    && BC.type.S==:Dirichlet) ? 1. : 0.
            NeuS   = (j==1    && BC.type.S==:Neumann  ) ? 1. : 0.
            inN    =  j==nc.y ? false  : true   
            DirN   = (j==nc.y && BC.type.N==:Dirichlet) ? 1. : 0.
            NeuN   = (j==nc.y && BC.type.N==:Neumann  ) ? 1. : 0.
            # Material coefficient
            kW = k.x[i,j]
            kE = k.x[i+1,j]
            kS = k.y[i,j]
            kN = k.y[i,j+1]
            ρ  = rho[i,j]
            Cp = cp[i,j]
            # Linear system coefficients
            if inS K[ii,iS] = kS .* (DirS + NeuS - 1) ./ dy .^ 2 end
            if inW K[ii,iW] = kW .* (DirW + NeuW - 1) ./ dx .^ 2 end
            K[ii,iC] = Cp .* ρ ./ Δt + (-kN .* (-DirN + NeuN - 1) ./ dy + kS .* (DirS - NeuS + 1) ./ dy) ./ dy + (-kE .* (-DirE + NeuE - 1) ./ dx + kW .* (DirW - NeuW + 1) ./ dx) ./ dx
            if inE K[ii,iE] = -kE .* (-DirE - NeuE + 1) ./ dx .^ 2 end
            if inN K[ii,iN] = -kN .* (-DirN - NeuN + 1) ./ dy .^ 2 end
        end
    end
    return flush!(K)
end

@doc raw"""
    BackwardEuler2Dc
"""
function BackwardEuler2Dc!(D, κ, Δx, Δy, Δt, ρ, cp, NC, BC, rhs, K, Num)
# dT/dt = kappa*d^2T_ij/dx_i^2 + Q_ij/rho/cp
# ----------------------------------------------------------------------- #
# Define coefficients ---
a   =   κ / Δx^2
b   =   κ / Δy^2
c   =   1 / Δt
# temp    =   copy(rhs)
# Multiply rhs with 1/Δt and add Q/ρ/cp ---    
# rhs  .= reshape(D.T,NC.x*NC.y).*c .+ reshape(D.Q,NC.x*NC.y)./ρ./cp
rhs  .= reshape(D.T,NC.x*NC.y).*c .+ 
            reshape(D.Q,NC.x*NC.y)./reshape(ρ,NC.x*NC.y)./cp

# Loop over the grid points ---
for i = 1:NC.x
    for j = 1:NC.y
        # Equation number ---
        ii          =   Num.T[i,j]
        # Stencil ---
        iS          =   ii - NC.x   # South
        iW          =   ii - 1      # West
        iC          =   ii          # Central
        iE          =   ii + 1      # East
        iN          =   ii + NC.x   # North
        # Boundaries ---
        # If an West index is required ---
        inW    =  i==1      ? false  : true
        DirW   = (i==1      && BC.type.W==:Dirichlet) ? 1. : 0.
        NeuW   = (i==1      && BC.type.W==:Neumann  ) ? 1. : 0.
        # If an East index is required ---
        inE    =  i==NC.x   ? false  : true
        DirE   = (i==NC.x   && BC.type.E==:Dirichlet) ? 1. : 0.
        NeuE   = (i==NC.x   && BC.type.E==:Neumann  ) ? 1. : 0.
        # If an South index is required
        inS    =  j==1      ? false  : true
        DirS   = (j==1      && BC.type.S==:Dirichlet) ? 1. : 0.
        NeuS   = (j==1      && BC.type.S==:Neumann  ) ? 1. : 0.
        # If an North index is required 
        inN    =  j==NC.y   ? false  : true
        DirN   = (j==NC.y   && BC.type.N==:Dirichlet) ? 1. : 0.
        NeuN   = (j==NC.y   && BC.type.N==:Neumann  ) ? 1. : 0.
        # Stencil ---
        if inS K[ii,iS]     = - b end
        if inW K[ii,iW]     = - a end
        K[ii,iC]            =   (2 + DirW + DirE - NeuW - NeuE)*a + (2 + DirS + DirN - NeuS - NeuN) *b + c
        if inE K[ii,iE]     = - a end    
        if inN K[ii,iN]     = - b end
        # Modify right hand side due to boundary conditions ------------- #
        rhs[ii]     +=  2*a*BC.val.W[j] * DirW +
                        2*a*BC.val.E[j] * DirE +
                        2*b*BC.val.S[i] * DirS +
                        2*b*BC.val.N[i] * DirN -
                        a*BC.val.W[j]*Δx * NeuW  + 
                        a*BC.val.E[j]*Δx * NeuE  - 
                        b*BC.val.S[i]*Δy * NeuS  + 
                        b*BC.val.N[i]*Δy * NeuN 
    end
end
# ------------------------------------------------------------------- #    
# Calculate temperature at new time step ---------------------------- #
D.T[:]  .=   K \ rhs[:]
D.T_ex[2:end-1,2:end-1]     .=    D.T
# ------------------------------------------------------------------- #
end

@doc raw"""
    CNA2Dc
"""
function CNA2Dc!(D, κ, Δx, Δy, Δt, ρ, cp, NC, BC, rhs, K1, K2, Num)
# dT/dt = kappa*d^2T_ij/dx_i^2 + Q_ij/rho/cp
# ----------------------------------------------------------------------- #

# Define coefficients ---
a       =   κ / 2 / Δx^2
b       =   κ / 2 / Δy^2
c       =   1 / Δt

# Loop over the grid points ---
for i = 1:NC.x
    for j = 1:NC.y
        # Equation number ---
        ii          =   Num.T[i,j]
        # Stencil ---
        iS          =   ii - NC.x   # South
        iW          =   ii - 1      # West
        iC          =   ii          # Central
        iE          =   ii + 1      # East
        iN          =   ii + NC.x   # North
        # Boundaries ---
        # If an West index is required ---
        inW    =  i==1    ? false  : true
        DirW   = (i==1    && BC.type.W==:Dirichlet) ? 1. : 0.
        NeuW   = (i==1    && BC.type.W==:Neumann  ) ? 1. : 0.
        # If an East index is required ---
        inE    =  i==NC.x ? false  : true
        DirE   = (i==NC.x && BC.type.E==:Dirichlet) ? 1. : 0.
        NeuE   = (i==NC.x && BC.type.E==:Neumann  ) ? 1. : 0.
        # If an South index is required
        inS    =  j==1      ? false  : true
        DirS   = (j==1      && BC.type.S==:Dirichlet) ? 1. : 0.
        NeuS   = (j==1      && BC.type.S==:Neumann  ) ? 1. : 0.
        # If an North index is required 
        inN    =  j==NC.y   ? false  : true
        DirN   = (j==NC.y   && BC.type.N==:Dirichlet) ? 1. : 0.
        NeuN   = (j==NC.y   && BC.type.N==:Neumann  ) ? 1. : 0.
        if inS
            K1[ii,iS]   =   - b
            K2[ii,iS]   =   b
        end
        if inE
            K1[ii,iE]   =   - a
            K2[ii,iE]   =   a            
        end
        K1[ii,iC]       =   c + (2 + DirW + DirE - NeuE - NeuW) * a + (2 + DirS + DirN - NeuS - NeuN) * b
        K2[ii,iC]       =   c - (2 + DirW + DirE - NeuE - NeuW) * a - (2 + DirS + DirN - NeuS - NeuN) * b
        if inW 
            K1[ii,iW]   =   - a
            K2[ii,iW]   =   a
        end           
        if inN
            K1[ii,iN]   =   - b
            K2[ii,iN]   =   b
        end
    end
end
# ------------------------------------------------------------------- #
# Berechnung der rechten Seite -------------------------------------- #
# rhs     .=   K2 * reshape(D.T,NC.x*NC.y) .+ reshape(D.Q,NC.x*NC.y)./ρ./cp
rhs     .=   K2 * reshape(D.T,NC.x*NC.y) .+ 
                reshape(D.Q,NC.x*NC.y)./reshape(ρ,NC.x*NC.y)./cp
# ------------------------------------------------------------------- #        
# Aenderung der rechten Seite durch die Randbedingungen ------------- #    
for i = 1:NC.x
    for j = 1:NC.y     
        ii      =   Num.T[i,j]
        # Boundaries         
        DirW    =   (i==1       && BC.type.W==:Dirichlet) ? 1. : 0.
        NeuW    =   (i==1       && BC.type.W==:Neumann  ) ? 1. : 0.
        DirE    =   (i==NC.x    && BC.type.E==:Dirichlet) ? 1. : 0.
        NeuE    =   (i==NC.x    && BC.type.E==:Neumann  ) ? 1. : 0.
        DirS    =   (j==1       && BC.type.S==:Dirichlet) ? 1. : 0.
        NeuS    =   (j==1       && BC.type.S==:Neumann  ) ? 1. : 0.
        DirN    =   (j==NC.y    && BC.type.N==:Dirichlet) ? 1. : 0.
        NeuN    =   (j==NC.y    && BC.type.N==:Neumann  ) ? 1. : 0.
        
        # Update right hand side ---
        rhs[ii]     = rhs[ii] +
                        4*a*BC.val.W[j] * DirW + 
                        4*a*BC.val.E[j] * DirE +
                        4*b*BC.val.S[i] * DirS +
                        4*b*BC.val.N[i] * DirN - 
                        2*a*BC.val.W[j]*Δx * NeuW +
                        2*a*BC.val.E[j]*Δx * NeuE -
                        2*b*BC.val.S[j]*Δy * NeuS + 
                        2*b*BC.val.N[j]*Δy * NeuN
    end
end
# ------------------------------------------------------------------- #    
# Compute new temperature ------------------------------------------- #
D.T[:]      .=    K1 \ rhs[:]
D.T_ex[2:end-1,2:end-1]     .=    D.T
# ------------------------------------------------------------------- #
end

@doc raw"""
    ADI2Dc
"""
function ADI2Dc!(T, κ, Δx, Δy, Δt, ρ, cp, NC, BC)
    # Function to solve 2D heat diffusion equation using the alternating direct
    # implicit finite difference scheme.
    # assuming constant k, rho, cp
    # dT/dt = kappa*d^2T_ij/dx_i^2 + Q_ij/rho/cp
    # ----------------------------------------------------------------------- #
    # Erstellung der durchlaufenden Indizes ----------------------------- #
    # Gleichungssystem fuer ADI Solver:
    Num     = (Th = reshape(1:NC.x*NC.y,NC.x,NC.y),
                Tv = reshape(1:NC.x*NC.y,NC.y,NC.x)')
    # Linear System of Equations -------------------------------------------- #
    ndof    =   maximum(Num.Th)
    A       =   ExtendableSparseMatrix(ndof,ndof)
    B       =   ExtendableSparseMatrix(ndof,ndof)
    C       =   ExtendableSparseMatrix(ndof,ndof)
    D       =   ExtendableSparseMatrix(ndof,ndof)
    rhs     =   zeros(ndof)
    temp    =   zeros(ndof)
    # Setup coefficient matrices ---------------------------------------- #
    a       =   κ / Δx^2
    b       =   κ / Δy^2
    c       =   1.0 / (Δt / 2.0)
    # First ADI step, vertical running scheme ---
    # Erster ADI Schritt: A*T^(l+1/2) = B*T^l -> vertical running scheme
    for i = 1:NC.x
        for j=1:NC.y
            # Equation number ---
            ii      =   Num.Tv[i,j]
            # Stencil ---
            iW      =   ii - NC.y
            iS      =   ii - 1
            iC      =   ii
            iN      =   ii + 1
            iE      =   ii + NC.y
            # Boundaries ---
            # If an West index is required ---
            inW    =  i==1    ? false  : true
            DirW   = (i==1    && BC.type.W==:Dirichlet) ? 1. : 0.
            NeuW   = (i==1    && BC.type.W==:Neumann  ) ? 1. : 0.
            # If an East index is required ---
            inE    =  i==NC.x ? false  : true
            DirE   = (i==NC.x && BC.type.E==:Dirichlet) ? 1. : 0.
            NeuE   = (i==NC.x && BC.type.E==:Neumann  ) ? 1. : 0.
            # If an South index is required
            inS    =  j==1      ? false  : true
            DirS   = (j==1      && BC.type.S==:Dirichlet) ? 1. : 0.
            NeuS   = (j==1      && BC.type.S==:Neumann  ) ? 1. : 0.
            # If an North index is required 
            inN    =  j==NC.y   ? false  : true
            DirN   = (j==NC.y   && BC.type.N==:Dirichlet) ? 1. : 0.
            NeuN   = (j==NC.y   && BC.type.N==:Neumann  ) ? 1. : 0.
            if inW
                B[ii,iW]   =  a 
            end
            if inN
                A[ii,iN]   =   - b
            end
            A[ii,iC]       =   c + (2 + DirS + DirN - NeuS - NeuN) * b
            B[ii,iC]       =   c - (2 + DirW + DirE - NeuW - NeuE) * a
            if inS
                A[ii,iS]   =   - b
            end
            if inE
                B[ii,iE]   =   a
            end        
        end
    end
    
    rhs  .=   B * reshape(T.T',(NC.y*NC.x,1)) .+ 
                    reshape(T.Q',(NC.y*NC.x,1))./reshape(ρ',(NC.y*NC.x,1))./cp

    # Update rhs from the boundary conditions ---
    for i = 1:NC.x
        for j = 1:NC.y
            ii  =   Num.Tv[i,j]
            # If an West index is required ---
            DirW   = (i==1    && BC.type.W==:Dirichlet) ? 1. : 0.
            NeuW   = (i==1    && BC.type.W==:Neumann  ) ? 1. : 0.
            # If an East index is required ---
            DirE   = (i==NC.x && BC.type.E==:Dirichlet) ? 1. : 0.
            NeuE   = (i==NC.x && BC.type.E==:Neumann  ) ? 1. : 0.
            # If an South index is required
            DirS   = (j==1      && BC.type.S==:Dirichlet) ? 1. : 0.
            NeuS   = (j==1      && BC.type.S==:Neumann  ) ? 1. : 0.
            # If an North index is required 
            DirN   = (j==NC.y   && BC.type.N==:Dirichlet) ? 1. : 0.
            NeuN   = (j==NC.y   && BC.type.N==:Neumann  ) ? 1. : 0.
            rhs[ii]     +=  2*a*BC.val.W[j] * DirW + 
                            2*a*BC.val.E[j] * DirE + 
                            2*b*BC.val.S[i] * DirS + 
                            2*b*BC.val.N[i] * DirN -
                            a*BC.val.W[j] * NeuW + 
                            a*BC.val.E[j] * NeuE - 
                            b*BC.val.S[i] * NeuS + 
                            b*BC.val.N[i] * NeuN    
        end
    end
    # Temperature at Δt/2 ---
    temp    .=  A \ rhs

    T.T     .=   reshape(temp,(NC.y,NC.x))'

    # Second ADI step, horizontal running scheme ---
    # Zweiter ADI Schritt: C*T^(l+1) = D*T^(l+1/2) -> horizontal running scheme
    for j = 1:NC.y
        for i=1:NC.x
            # Equation number ---
            ii      =   Num.Th[i,j]
            # Stencil ---
            iS      =   ii - NC.x
            iW      =   ii - 1        
            iC      =   ii
            iE      =   ii + 1
            iN      =   ii + NC.x
            # Boundaries ---
            # If an West index is required ---
            inW    =  i==1    ? false  : true
            DirW   = (i==1    && BC.type.W==:Dirichlet) ? 1. : 0.
            NeuW   = (i==1    && BC.type.W==:Neumann  ) ? 1. : 0.
            # If an East index is required ---
            inE    =  i==NC.x ? false  : true
            DirE   = (i==NC.x && BC.type.E==:Dirichlet) ? 1. : 0.
            NeuE   = (i==NC.x && BC.type.E==:Neumann  ) ? 1. : 0.
            # If an South index is required
            inS    =  j==1      ? false  : true
            DirS   = (j==1      && BC.type.S==:Dirichlet) ? 1. : 0.
            NeuS   = (j==1      && BC.type.S==:Neumann  ) ? 1. : 0.
            # If an North index is required 
            inN    =  j==NC.y   ? false  : true
            DirN   = (j==NC.y   && BC.type.N==:Dirichlet) ? 1. : 0.
            NeuN   = (j==NC.y   && BC.type.N==:Neumann  ) ? 1. : 0.
            if inN
                D[ii,iN]   =  b
            end
            if inW
                C[ii,iW]   =  - a
            end
            C[ii,iC]       =   c + (2 + DirW + DirE - NeuW - NeuE) * a
            D[ii,iC]       =   c - (2 + DirS + DirN - NeuS - NeuN) * b
            if inE
                C[ii,iE]   =   - a
            end   
            if inS
                D[ii,iS]   =   b
            end        
        end
    end

    # Update rhs to T^{n+1/2} --- 
    rhs  .=   D * reshape(T.T,(NC.y*NC.x,1)) .+ 
                    reshape(T.Q,(NC.y*NC.x,1))./reshape(ρ,(NC.y*NC.x,1))./cp
    # rhs  .=   D * temp .+ 
    #                 reshape(T.Q,(NC.y*NC.x,1))./reshape(ρ,(NC.y*NC.x,1))./cp
    
    # Update rhs from the boundary conditions ---
    for j = 1:NC.y
        for i = 1:NC.x
            ii  =   Num.Th[i,j]
            # If an West index is required ---
            DirW   = (i==1    && BC.type.W==:Dirichlet) ? 1. : 0.
            NeuW   = (i==1    && BC.type.W==:Neumann  ) ? 1. : 0.
            # If an East index is required ---
            DirE   = (i==NC.x && BC.type.E==:Dirichlet) ? 1. : 0.
            NeuE   = (i==NC.x && BC.type.E==:Neumann  ) ? 1. : 0.
            # If an South index is required
            DirS   = (j==1      && BC.type.S==:Dirichlet) ? 1. : 0.
            NeuS   = (j==1      && BC.type.S==:Neumann  ) ? 1. : 0.
            # If an North index is required 
            DirN   = (j==NC.y   && BC.type.N==:Dirichlet) ? 1. : 0.
            NeuN   = (j==NC.y   && BC.type.N==:Neumann  ) ? 1. : 0.
            rhs[ii]     +=  2*a*BC.val.W[j] * DirW + 
                            2*a*BC.val.E[j] * DirE + 
                            2*b*BC.val.S[i] * DirS + 
                            2*b*BC.val.N[i] * DirN -
                            a*BC.val.W[j] * NeuW + 
                            a*BC.val.E[j] * NeuE - 
                            b*BC.val.S[i] * NeuS + 
                            b*BC.val.N[i] * NeuN    
        end
    end

    # Temperature at Δt ---
    temp    .=  C \ rhs

    T.T     .=   reshape(temp,(NC.x,NC.y))

    T.T_ex[2:end-1,2:end-1]     .=    T.T
end

@doc raw"""
    Poisson2Dc
"""
function Poisson2Dc!(D,NC,P,BC,Δ,K,rhs,Num)
# Function to solve 2D heat diffusion equation using the explicit finite
# difference scheme
# [Q] = W/m^3
# ----------------------------------------------------------------------- #
    
    # a       =   1.0 / Δ.x[1]^2.0
    # b       =   1.0 / Δ.y[1]^2.0

    #  --------------------------------------------- #
    rhs     .=   - reshape(D.Q, NC.x*NC.y, 1) ./ P.k

    for i=1:NC.x
        for j=1:NC.y
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
    end

    D.T[:]  .=   K \ rhs[:]
    
end

@doc raw"""
    Poisson2D
"""
function Poisson2D!( T, Q, kx, ky, Δx, Δy, NC, BC, K, rhs, Num ) 
    #  --------------------------------------------- #
    rhs     .=   - reshape(Q, NC.x*NC.y, 1)

    for i=1:NC.x
        for j=1:NC.y
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
    end

    T[:]  .=   K \ rhs[:]
end
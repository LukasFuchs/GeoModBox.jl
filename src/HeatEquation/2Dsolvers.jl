using ExtendableSparse

# ======================================================================= #
# Time-dependent solvers, constant thermal parameters =================== #
# ======================================================================= #
"""
    ForwardEuler2Dc!(D, κ, Δx, Δy, Δt, ρ, cp, NC, BC)

Solves the two dimensional heat diffusion equation assuming constant thermal 
parameters using an explicit, forward Euler finite difference scheme.

The temperature is defined on central nodes and the heat flux on the vertices. 
Boundary conditions are currently limited to Dirichlet and Neumann. Using central 
temperature nodes requires external ghost nodes, which are used to define the 
boundary conditions. 

    D           : Tuple, containing the regular temperature array T and 
                  the extended array containing the ghost nodes T_ex
    κ           : Diffusivity [ m²/s ]
    Δx          : Horizontal grid spacing [ m ]
    Δy          : Vertical grid spacing [ m ]
    Δt          : Time step [ s ]
    ρ           : Density [ kg/m³ ]
    cp          : Specific heat capacity [ J/kg/K ]
    NC          : Tuple containing the numer of centroids in x- and y-direction
    BC          : Tuple for the boundary condition

Optional input values (to include a heat source): 
    Q           : Volumetric heat production rate [ W/m^3 ]
    ρ₀          : Reference density [ kg/m^3 ]
    cp          : Specific heat capacity [ J/kg/K ]

"""
function ForwardEuler2Dc!(D, κ, Δx, Δy, Δt, NC, BC; 
                Q = zeros(NC...), ρ₀ = 3300.0, cp = 1200.0 )
    # Function to solve 2D heat diffusion equation using the explicit finite
    # difference scheme
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
                Q[i,j] * Δt / ρ₀ / cp
        end
    end
    # ------------------------------------------------------------------- #
    # Update extended temperature --------------------------------------- #
    D.T_ex[2:end-1,2:end-1]     .=  D.T
    # ------------------------------------------------------------------- #    
end

"""
    ComputeResiduals2Dc!(R, T, T_ex, T0, T_ex0, ∂2T, Q, κ, BC, Δ, Δt;C=0)

Function to calculate the residual of the two dimensional heat diffusion equation assuming 
constant thermal parameters and radiogenic heating only. The residual is required for 
defection correction to solve the linear system of equations. 

The temperature is defined on central nodes and the heat flux on the vertices. 
Boundary conditions are currently limited to Dirichlet and Neumann. Using central 
temperature nodes requires external ghost nodes, which are used to define the 
boundary conditions. The caluclated residual is used in an iteration loop to calculate 
the correction term of the initial temperature guess. The coefficient matrix is build 
via `AssembleMatrix2Dc()`. 

    R           : 2D array for the residual
    T           : 2D centroid temperature field of the next time step
    T_ex        : 2D extended centroid temperature field including the ghost nodes
    T_0         : 2D centroid temperature field of the current time step
    T_ex0       : 2D extended centroid temperature field of the current time step
    ∂2T         : Tuple containing the second space derivatives of the temperature
                  in the horizontal and vertical direction
    Q           : Volumetric heat production rate [ W/m³ ]
    κ           : Diffusivity [ m²/s ]
    BC          : Tuple for the boundary condition
    Δ           : Tuple or strucutre containint the horizontal and vertical grid resolution [ m ]
    Δt          : Time step [ s ]

Optional input values: 
    C           : Constant defining the residual for a certain finite difference discretization 
                  method, i.e.: 
                        C = 0   -> implicit, backward Euler discretization (default)
                        C = 0.5 -> Crank-Nicolson discretization
                        C = 1   -> explicit, forward Euler discretization
    Q           : Volumetric heat production rate [ W/m^3 ]
    ρ₀          : Reference density [ kg/m^3 ]
    cp          : Specific heat capacity [ J/kg/K ]
    
"""
function ComputeResiduals2Dc!( R, T, T_ex, T0, T_ex0, ∂2T, κ, BC, Δ, Δt;
                C = 0, Q = 0.0, ρ₀ = 3300.0, cp = 1200.0 )
    if C < 1
        # Implicit 
        @. T_ex[2:end-1,2:end-1] = T 
        @. T_ex[  1,2:end-1] = (BC.type.W==:Dirichlet) * (2*BC.val.W - T_ex[    2,2:end-1]) + (BC.type.W==:Neumann) * (T_ex[    2,2:end-1] - Δ.x*BC.val.W)
        @. T_ex[end,2:end-1] = (BC.type.E==:Dirichlet) * (2*BC.val.E - T_ex[end-1,2:end-1]) + (BC.type.E==:Neumann) * (T_ex[end-1,2:end-1] + Δ.x*BC.val.E)
        @. T_ex[2:end-1,  1] = (BC.type.S==:Dirichlet) * (2*BC.val.S - T_ex[2:end-1,    2]) + (BC.type.S==:Neumann) * (T_ex[2:end-1,    2] - Δ.y*BC.val.S)
        @. T_ex[2:end-1,end] = (BC.type.N==:Dirichlet) * (2*BC.val.N - T_ex[2:end-1,end-1]) + (BC.type.N==:Neumann) * (T_ex[2:end-1,end-1] + Δ.y*BC.val.N)
        @. ∂2T.∂x2  = (T_ex[1:end-2,2:end-1] - 2.0 * T_ex[2:end-1,2:end-1] + T_ex[3:end,2:end-1])/Δ.x/Δ.x
        @. ∂2T.∂y2  = (T_ex[2:end-1,1:end-2] - 2.0 * T_ex[2:end-1,2:end-1] + T_ex[2:end-1,3:end])/Δ.y/Δ.y
        if C==0.5
            # CNA
            @. T_ex0[2:end-1,2:end-1] = T0 
            @. T_ex0[  1,2:end-1] = (BC.type.W==:Dirichlet) * (2*BC.val.W - T_ex0[    2,2:end-1]) + (BC.type.W==:Neumann) * (T_ex0[    2,2:end-1] - Δ.x*BC.val.W)
            @. T_ex0[end,2:end-1] = (BC.type.E==:Dirichlet) * (2*BC.val.E - T_ex0[end-1,2:end-1]) + (BC.type.E==:Neumann) * (T_ex0[end-1,2:end-1] + Δ.x*BC.val.E)
            @. T_ex0[2:end-1,  1] = (BC.type.S==:Dirichlet) * (2*BC.val.S - T_ex0[2:end-1,    2]) + (BC.type.S==:Neumann) * (T_ex0[2:end-1,    2] - Δ.y*BC.val.S)
            @. T_ex0[2:end-1,end] = (BC.type.N==:Dirichlet) * (2*BC.val.N - T_ex0[2:end-1,end-1]) + (BC.type.N==:Neumann) * (T_ex0[2:end-1,end-1] + Δ.y*BC.val.N)
            @. ∂2T.∂x20  = (T_ex0[1:end-2,2:end-1] - 2.0 * T_ex0[2:end-1,2:end-1] + T_ex0[3:end,2:end-1])/Δ.x/Δ.x
            @. ∂2T.∂y20  = (T_ex0[2:end-1,1:end-2] - 2.0 * T_ex0[2:end-1,2:end-1] + T_ex0[2:end-1,3:end])/Δ.y/Δ.y
        end
    else
        # explicit
        @. T_ex0[2:end-1,2:end-1] = T0 
        @. T_ex0[  1,2:end-1] = (BC.type.W==:Dirichlet) * (2*BC.val.W - T_ex0[    2,2:end-1]) + (BC.type.W==:Neumann) * (T_ex0[    2,2:end-1] - Δ.x*BC.val.W)
        @. T_ex0[end,2:end-1] = (BC.type.E==:Dirichlet) * (2*BC.val.E - T_ex0[end-1,2:end-1]) + (BC.type.E==:Neumann) * (T_ex0[end-1,2:end-1] + Δ.x*BC.val.E)
        @. T_ex0[2:end-1,  1] = (BC.type.S==:Dirichlet) * (2*BC.val.S - T_ex0[2:end-1,    2]) + (BC.type.S==:Neumann) * (T_ex0[2:end-1,    2] - Δ.y*BC.val.S)
        @. T_ex0[2:end-1,end] = (BC.type.N==:Dirichlet) * (2*BC.val.N - T_ex0[2:end-1,end-1]) + (BC.type.N==:Neumann) * (T_ex0[2:end-1,end-1] + Δ.y*BC.val.N)
        @. ∂2T.∂x20  = (T_ex0[1:end-2,2:end-1] - 2.0 * T_ex0[2:end-1,2:end-1] + T_ex0[3:end,2:end-1])/Δ.x/Δ.x
        @. ∂2T.∂y20  = (T_ex0[2:end-1,1:end-2] - 2.0 * T_ex0[2:end-1,2:end-1] + T_ex0[2:end-1,3:end])/Δ.y/Δ.y
    end
    @. R     = (T - T0)/Δt - κ*((1-C)*(∂2T.∂x2 + ∂2T.∂y2) + C*(∂2T.∂x20 + ∂2T.∂y20)) - Q/ρ₀/cp
end

"""
    AssembleMatrix2Dc( κ, BC, Num, nc, Δ, Δt )

Function to build the coefficient matrix K for the unknown centroid temperature field in the 
2D heat diffusion equation assuming constant thermal parameter and radiogenig heating only. 

The coefficient matrix is build using a five-point finite difference stencil, resulting in 
a five, non-zero diagonal matrix to solve the system of equations. 

    κ       : Diffusivity [ m²/s ]
    BC      : Tuple for the boundary condition
    Num     : Tuple or structure containing the global numbering of the centroids
    nc      : Tuple or structure containing the number of centroids in the horizontal and vertical direction 
    Δ       : Tuple or structure containint the horizontal and vertical grid resolution
    Δt      : Time step [ s ]
"""
function AssembleMatrix2Dc( κ, BC, Num, nc, Δ, Δt;C=0 )
    # Linear system of equation
    ndof   = maximum(Num.T)
    K      = ExtendableSparseMatrix(ndof, ndof)
    # dx, dy = Δ.x, Δ.y
    # Define coefficients ---
    a   =   (κ*(1-C)) / Δ.x^2
    b   =   (κ*(1-C)) / Δ.y^2
    c   =   1 / Δt
    #############################
    #       Heat equation       #
    #############################
    for i = 1:nc.x
        for j = 1:nc.y
            # Equation number ---
            ii          =   Num.T[i,j]
            # Stencil ---
            iS          =   ii - nc.x   # South
            iW          =   ii - 1      # West
            iC          =   ii          # Central
            iE          =   ii + 1      # East
            iN          =   ii + nc.x   # North
            # Boundaries ---
            # If an West index is required ---
            inW    =  i==1      ? false  : true
            DirW   = (i==1      && BC.type.W==:Dirichlet) ? 1. : 0.
            NeuW   = (i==1      && BC.type.W==:Neumann  ) ? 1. : 0.
            # If an East index is required ---
            inE    =  i==nc.x   ? false  : true
            DirE   = (i==nc.x   && BC.type.E==:Dirichlet) ? 1. : 0.
            NeuE   = (i==nc.x   && BC.type.E==:Neumann  ) ? 1. : 0.
            # If an South index is required
            inS    =  j==1      ? false  : true
            DirS   = (j==1      && BC.type.S==:Dirichlet) ? 1. : 0.
            NeuS   = (j==1      && BC.type.S==:Neumann  ) ? 1. : 0.
            # If an North index is required 
            inN    =  j==nc.y   ? false  : true
            DirN   = (j==nc.y   && BC.type.N==:Dirichlet) ? 1. : 0.
            NeuN   = (j==nc.y   && BC.type.N==:Neumann  ) ? 1. : 0.
            # Stencil ---
            if inS K[ii,iS]     = - b end
            if inW K[ii,iW]     = - a end
            K[ii,iC]            =   (2 + DirW + DirE - NeuW - NeuE)*a + (2 + DirS + DirN - NeuS - NeuN) *b + c
            if inE K[ii,iE]     = - a end    
            if inN K[ii,iN]     = - b end
        end
    end
    return flush!(K)
end

"""
    BackwardEuler2Dc!(D, κ, Δx, Δy, Δt, ρ, cp, NC, BC, rhs, K, Num)

Solves the two dimensional heat diffusion equation assuming constant thermal 
parameters using an implicit, backward Euler finite difference scheme for a 
linear problem, i.e. a left-matrix division. 

The temperature is defined on central nodes and the heat flux on the vertices. 
Boundary conditions are currently limited to Dirichlet and Neumann. Boundary conditions
are directly implemented within the system of equations.  

    D           : Tuple, containing the regular temperature array T and 
                  the extended array containing the ghost nodes T_ex
    κ           : Diffusivity [ m²/s ]
    Δx          : Horizontal grid spacing [ m ]
    Δy          : Vertical grid spacing [ m ]
    Δt          : Time step [ s ]
    ρ           : Density [ kg/m³ ]
    cp          : Specific heat capacity [ J/kg/K ]
    NC          : Tuple containing the numer of centroids in x- and y-direction
    BC          : Tuple for the boundary condition
    rhs         : Known right-hand side vector
    K           : Coefficient matrix in sparse format
    Num         : Tuple or structure containing the global centroid numbering

Optional input values (to include a heat source): 
    Q           : Volumetric heat production rate [ W/m^3 ]
    ρ₀          : Reference density [ kg/m^3 ]
    cp          : Specific heat capacity [ J/kg/K ]

"""
function BackwardEuler2Dc!(D, κ, Δx, Δy, Δt, NC, BC, rhs, K, Num; 
                        Q = zeros(NC...), ρ₀ = 3300.0, cp = 1200.0 )
# dT/dt = kappa*d^2T_ij/dx_i^2 + Q_ij/ρ/cp
# ----------------------------------------------------------------------- #
# Define coefficients ---
a   =   κ / Δx^2
b   =   κ / Δy^2
c   =   1 / Δt

rhs  .= reshape(D.T,NC.x*NC.y).*c .+ 
            reshape(Q,NC.x*NC.y)./ρ₀./cp

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

"""
    CNA2Dc!(D, κ, Δx, Δy, Δt, ρ, cp, NC, BC, rhs, K1, K2, Num)

Solves the two dimensional heat diffusion equation assuming constant thermal 
parameters using the Crank-Nicolson finite difference scheme for a 
linear problem, i.e. a left-matrix division. 

The temperature is defined on central nodes and the heat flux on the vertices. 
Boundary conditions are currently limited to Dirichlet and Neumann. Boundary conditions
are directly implemented within the system of equations.  

    D           : Tuple, containing the regular temperature array T and 
                  the extended array containing the ghost nodes T_ex
    κ           : Diffusivity [ m²/s ]
    Δx          : Horizontal grid spacing [ m ]
    Δy          : Vertical grid spacing [ m ]
    Δt          : Time step [ s ]
    ρ           : Density [ kg/m³ ]
    cp          : Specific heat capacity [ J/kg/K ]
    NC          : Tuple containing the numer of centroids in x- and y-direction
    BC          : Tuple for the boundary condition
    rhs         : Known right-hand side vector
    K1          : Coefficient matrix in sparse format for the unknown temperature
    K1          : Coefficient matrix in sparse format for the known temperature
    Num         : Tuple or structure containing the global centroid numbering

Optional input values (to include a heat source): 
    Q           : Volumetric heat production rate [ W/m^3 ]
    ρ₀          : Reference density [ kg/m^3 ]
    cp          : Specific heat capacity [ J/kg/K ]
"""
function CNA2Dc!(D, κ, Δx, Δy, Δt, NC, BC, rhs, K1, K2, Num; 
                Q = zeros(NC...), ρ₀ = 3300.0, cp = 1200.0 )
# dT/dt = kappa*d^2T_ij/dx_i^2 + Q_ij/ρ/cp
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
rhs     .=   K2 * reshape(D.T,NC.x*NC.y) .+ 
                reshape(Q,NC.x*NC.y)./ρ₀./cp
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

"""
    ADI2Dc!(T, κ, Δx, Δy, Δt, ρ, cp, NC, BC)

Solves the two dimensional heat diffusion equation assuming constant thermal 
parameters using the alternating-direction implicit finite difference scheme for a 
linear problem, i.e. a left-matrix division. 

The temperature is defined on central nodes and the heat flux on the vertices. 
Boundary conditions are currently limited to Dirichlet and Neumann. Boundary conditions
are directly implemented within the system of equations.  

    T           : Tuple, containing the regular temperature array T, 
                  the extended array containing the ghost nodes T_ex, and 
                  the volumentric heat production rate Q 
    κ           : Diffusivity [ m²/s ]
    Δx          : Horizontal grid spacing [ m ]
    Δy          : Vertical grid spacing [ m ]
    Δt          : Time step [ s ]
    ρ           : Density [ kg/m³ ]
    cp          : Specific heat capacity [ J/kg/K ]
    NC          : Tuple containing the numer of centroids in x- and y-direction
    BC          : Tuple for the boundary condition

Optional input values (to include a heat source): 
    Q           : Volumetric heat production rate [ W/m^3 ]
    ρ₀          : Reference density [ kg/m^3 ]
    cp          : Specific heat capacity [ J/kg/K ]
"""
function ADI2Dc!(T, κ, Δx, Δy, Δt, NC, BC; 
                Q = zeros(NC...), ρ₀ = 3300.0, cp = 1200.0 )
    # Function to solve 2D heat diffusion equation using the alternating direct
    # implicit finite difference scheme.
    # assuming constant k, ρ, cp
    # dT/dt = kappa*d^2T_ij/dx_i^2 + Q_ij/ρ/cp
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
                    reshape(Q',(NC.y*NC.x,1))./ρ₀./cp

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
                    reshape(Q,(NC.y*NC.x,1))./ρ₀./cp
    
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

# ======================================================================= #
# Time-dependent solvers, variable thermal parameters =================== #
# ======================================================================= #
"""
    ComputeResiduals2D!(R, T, T_ex, T0, T_ex0, Q, ∂T, q, ρ, Cp, k, BC, Δ, Δt;C=0)

Function to calculate the residual of the two dimensional heat diffusion equation assuming 
variable thermal parameters and radiogenic heating only. The residual is required for 
defection correction to solve the linear system of equations. 

The temperature is defined on central nodes and the heat flux on the vertices. 
Boundary conditions are currently limited to Dirichlet and Neumann. Using central 
temperature nodes requires external ghost nodes, which are used to define the 
boundary conditions. The caluclated residual is used in an iteration loop to calculate 
the correction term of the initial temperature guess. The coefficient matrix is build 
via `AssembleMatrix2D()`. 

    R           : 2D array for the residual
    T           : 2D centroid temperature field of the next time step
    T_ex        : 2D extended centroid temperature field including the ghost nodes
    T_0         : 2D centroid temperature field of the current time step
    T_ex0       : 2D extended centroid temperature field of the current time step
    Q           : Volumetric heat production rate [ W/m³ ]
    ∂T          : Tuple containing the first space derivatives of the temperature
                  in the horizontal and vertical direction
    q           : Tuple or structure containing the 2D horizontal and vertical heat flux
    ρ           : Density [ kg/m³ ]
    Cp          : Specific heat capacity [ J/kg/K ]
    k           : Thermal conductivity [ W/m/K ]
    BC          : Tuple for the boundary condition
    Δ           : Tuple or strucutre containint the horizontal and vertical grid resolution [ m ]
    Δt          : Time step [ s ]

Optional input values: 
    C           : Constant defining the residual for a certain finite difference discretization 
                  method, i.e.: 
                        C = 0   -> implicit, backward Euler discretization (default)
                        C = 0.5 -> Crank-Nicolson discretization
                        C = 1   -> explicit, forward Euler discretization
    
"""
function ComputeResiduals2D!(R, T, T_ex, T0, T_ex0, Q, ∂T, q, ρ, Cp, k, BC, Δ, Δt;C=0)
    if C < 1
        @. T_ex[2:end-1,2:end-1] = T 
        @. T_ex[  1,2:end-1] = (BC.type.W==:Dirichlet) * (2*BC.val.W - T_ex[    2,2:end-1]) + (BC.type.W==:Neumann) * (T_ex[    2,2:end-1] - Δ.x/k.x[  1,:]*BC.val.W)
        @. T_ex[end,2:end-1] = (BC.type.E==:Dirichlet) * (2*BC.val.E - T_ex[end-1,2:end-1]) + (BC.type.E==:Neumann) * (T_ex[end-1,2:end-1] + Δ.x/k.x[end,:]*BC.val.E)
        @. T_ex[2:end-1,  1] = (BC.type.S==:Dirichlet) * (2*BC.val.S - T_ex[2:end-1,    2]) + (BC.type.S==:Neumann) * (T_ex[2:end-1,    2] - Δ.y/k.y[:,  1]*BC.val.S)
        @. T_ex[2:end-1,end] = (BC.type.N==:Dirichlet) * (2*BC.val.N - T_ex[2:end-1,end-1]) + (BC.type.N==:Neumann) * (T_ex[2:end-1,end-1] + Δ.y/k.y[:,end]*BC.val.N)
        @. ∂T.∂x = (T_ex[2:end,2:end-1] - T_ex[1:end-1,2:end-1])/Δ.x
        @. ∂T.∂y = (T_ex[2:end-1,2:end] - T_ex[2:end-1,1:end-1])/Δ.y
        @. q.x   = -k.x * ∂T.∂x
        @. q.y   = -k.y * ∂T.∂y
        if C==0.5
            @. T_ex0[2:end-1,2:end-1] = T0 
            @. T_ex0[  1,2:end-1] = (BC.type.W==:Dirichlet) * (2*BC.val.W - T_ex0[    2,2:end-1]) + (BC.type.W==:Neumann) * (T_ex0[    2,2:end-1] - Δ.x/k.x[  1,:]*BC.val.W)
            @. T_ex0[end,2:end-1] = (BC.type.E==:Dirichlet) * (2*BC.val.E - T_ex0[end-1,2:end-1]) + (BC.type.E==:Neumann) * (T_ex0[end-1,2:end-1] + Δ.x/k.x[end,:]*BC.val.E)
            @. T_ex0[2:end-1,  1] = (BC.type.S==:Dirichlet) * (2*BC.val.S - T_ex0[2:end-1,    2]) + (BC.type.S==:Neumann) * (T_ex0[2:end-1,    2] - Δ.y/k.y[:,  1]*BC.val.S)
            @. T_ex0[2:end-1,end] = (BC.type.N==:Dirichlet) * (2*BC.val.N - T_ex0[2:end-1,end-1]) + (BC.type.N==:Neumann) * (T_ex0[2:end-1,end-1] + Δ.y/k.y[:,end]*BC.val.N)
            @. ∂T.∂x = (T_ex0[2:end,2:end-1] - T_ex0[1:end-1,2:end-1])/Δ.x
            @. ∂T.∂y = (T_ex0[2:end-1,2:end] - T_ex0[2:end-1,1:end-1])/Δ.y
            @. q.x0  = -k.x * ∂T.∂x
            @. q.y0  = -k.y * ∂T.∂y
        end
    else
        @. T_ex0[2:end-1,2:end-1] = T0 
        @. T_ex0[  1,2:end-1] = (BC.type.W==:Dirichlet) * (2*BC.val.W - T_ex0[    2,2:end-1]) + (BC.type.W==:Neumann) * (T_ex0[    2,2:end-1] - Δ.x/k.x[  1,:]*BC.val.W)
        @. T_ex0[end,2:end-1] = (BC.type.E==:Dirichlet) * (2*BC.val.E - T_ex0[end-1,2:end-1]) + (BC.type.E==:Neumann) * (T_ex0[end-1,2:end-1] + Δ.x/k.x[end,:]*BC.val.E)
        @. T_ex0[2:end-1,  1] = (BC.type.S==:Dirichlet) * (2*BC.val.S - T_ex0[2:end-1,    2]) + (BC.type.S==:Neumann) * (T_ex0[2:end-1,    2] - Δ.y/k.y[:,  1]*BC.val.S)
        @. T_ex0[2:end-1,end] = (BC.type.N==:Dirichlet) * (2*BC.val.N - T_ex0[2:end-1,end-1]) + (BC.type.N==:Neumann) * (T_ex0[2:end-1,end-1] + Δ.y/k.y[:,end]*BC.val.N)
        @. ∂T.∂x = (T_ex0[2:end,2:end-1] - T_ex0[1:end-1,2:end-1])/Δ.x
        @. ∂T.∂y = (T_ex0[2:end-1,2:end] - T_ex0[2:end-1,1:end-1])/Δ.y
        @. q.x0  = -k.x * ∂T.∂x
        @. q.y0  = -k.y * ∂T.∂y
    end
    # @. R     = (T - T0)/Δt - ((1-C)*(∂2T.∂x2 + ∂2T.∂y2) + C*(∂2T.∂x20 + ∂2T.∂y20)) - Q
    @. R     = ρ*Cp*(T - T0)/Δt + 
                    (1-C)*((q.x[2:end,:] - q.x[1:end-1,:])/Δ.x + (q.y[:,2:end] - q.y[:,1:end-1])/Δ.y) + 
                    C*((q.x0[2:end,:] - q.x0[1:end-1,:])/Δ.x + (q.y0[:,2:end] - q.y0[:,1:end-1])/Δ.y) - 
                    Q

    # @. T_ex[2:end-1,2:end-1] = T 
    # @. T_ex[  1,2:end-1] = (BC.type.W==:Dirichlet) * (2*BC.val.W - T_ex[    2,2:end-1]) + (BC.type.W==:Neumann) * (T_ex[    2,2:end-1] - Δ.x/k.x[  1,:]*BC.val.W)
    # @. T_ex[end,2:end-1] = (BC.type.E==:Dirichlet) * (2*BC.val.E - T_ex[end-1,2:end-1]) + (BC.type.E==:Neumann) * (T_ex[end-1,2:end-1] + Δ.x/k.x[end,:]*BC.val.E)
    # @. T_ex[2:end-1,  1] = (BC.type.S==:Dirichlet) * (2*BC.val.S - T_ex[2:end-1,    2]) + (BC.type.S==:Neumann) * (T_ex[2:end-1,    2] - Δ.y/k.y[:,  1]*BC.val.S)
    # @. T_ex[2:end-1,end] = (BC.type.N==:Dirichlet) * (2*BC.val.N - T_ex[2:end-1,end-1]) + (BC.type.N==:Neumann) * (T_ex[2:end-1,end-1] + Δ.y/k.y[:,end]*BC.val.N)
    # @. ∂T.∂x = (T_ex[2:end,2:end-1] - T_ex[1:end-1,2:end-1])/Δ.x
    # @. ∂T.∂y = (T_ex[2:end-1,2:end] - T_ex[2:end-1,1:end-1])/Δ.y
    # @. q.x   = -k.x * ∂T.∂x
    # @. q.y   = -k.y * ∂T.∂y
    # @. R     = ρ*Cp*(T - T0)/Δt + (q.x[2:end,:] - q.x[1:end-1,:])/Δ.x + (q.y[:,2:end] - q.y[:,1:end-1])/Δ.y - Q
end

"""
    AssembleMatrix2D(ρ, cp, k, BC, Num, nc, Δ, Δt;C=0)

Function to build the coefficient matrix K for the unknown centroid temperature field in the 
2D heat diffusion equation assuming variable thermal parameter and radiogenig heating only. 

The coefficient matrix is build using a five-point finite difference stencil, resulting in 
a five, non-zero diagonal matrix to solve the system of equations. 

    ρ       : Density [ kg/m³ ]
    cp      : Specific heat capacity [ J/kg/K ]
    k       : Thermal conductivity [ W/m/K ]
    BC      : Tuple for the boundary condition
    Num     : Tuple or structure containing the global numbering of the centroids
    nc      : Tuple or structure containing the number of centroids in the horizontal and vertical direction 
    Δ       : Tuple or structure containint the horizontal and vertical grid resolution
    Δt      : Time step [ s ]

Optional input values: 
    C           : Constant defining the residual for a certain finite difference discretization 
                  method, i.e.: 
                        C = 0   -> implicit, backward Euler discretization (default)
                        C = 0.5 -> Crank-Nicolson discretization
                        C = 1   -> explicit, forward Euler discretization
"""
function AssembleMatrix2D(ρ, cp, k, BC, Num, nc, Δ, Δt;C=0)
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
            kW = k.x[i,j]*(1-C)
            kE = k.x[i+1,j]*(1-C)
            kS = k.y[i,j]*(1-C)
            kN = k.y[i,j+1]*(1-C)
            # ρ  = ρ[i,j]
            # Cp = cp[i,j]
            # Linear system coefficients
            if inS K[ii,iS] = kS .* (DirS + NeuS - 1) ./ dy .^ 2 end
            if inW K[ii,iW] = kW .* (DirW + NeuW - 1) ./ dx .^ 2 end
            K[ii,iC] = cp[i,j] .* ρ[i,j] ./ Δt + (-kN .* (-DirN + NeuN - 1) ./ dy + kS .* (DirS - NeuS + 1) ./ dy) ./ dy + (-kE .* (-DirE + NeuE - 1) ./ dx + kW .* (DirW - NeuW + 1) ./ dx) ./ dx
            if inE K[ii,iE] = -kE .* (-DirE - NeuE + 1) ./ dx .^ 2 end
            if inN K[ii,iN] = -kN .* (-DirN - NeuN + 1) ./ dy .^ 2 end
        end
    end
    return flush!(K)
end

# ======================================================================= #
# Steady state solvers, constant thermal parameters ===================== #
# ======================================================================= #
"""
    Poisson2Dc!(D,NC,P,BC,Δ,K,rhs,Num)

Function to solve the two dimensional poisson equation of the steady state 
heat diffusion equation assuming constant thermal parameters and radiogenic 
heat sources only. 

The temperature is defined on central nodes and the heat flux on the vertices. 
Boundary conditions are currently limited to Dirichlet and Neumann. Boundary conditions
are directly implemented within the system of equations.  

    D       : Tuple or structure containing the centroid temperature field 
              and volumetic heat production rate. 
    NC      : Tuple or structure containing the number of horizontal and vertical 
              centroids.
    P       : Tuple or structure containing the thermal conductivity [ W/m/K ].
    BC      : Tuple or structure containint the boundary conditions.
    Δ       : Tuple or strucutre containint the horizontal and vertical grid resolution [ m ]
    K       : Coefficient matrix in sparse format
    rhs     : Known right-hand side vector
    Num     : Tuple or structure containing the global centroid numbering
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

# ======================================================================= #
# Steady state solvers, variable thermal parameters ===================== #
# ======================================================================= #
"""
    Poisson2D!( T, Q, kx, ky, Δx, Δy, NC, BC, K, rhs, Num ) 

Function to solve the two dimensional poisson equation of the steady state 
heat diffusion equation assuming variable thermal parameters and radiogenic 
heat sources only. 

The temperature is defined on central nodes and the heat flux on the vertices. 
Boundary conditions are currently limited to Dirichlet and Neumann. Boundary conditions
are directly implemented within the system of equations.  

    T       : Tuple or structure containing the centroid temperature field. 
    Q       : Tuple or structure containing the centroid volumentric heat prodcution rate [ W/m³ ].
    kx      : Thermal conductivity in the horizontal direction [ W/m/K ].
    ky      : Thermal conductivity in the vertical direction [ W/m/K ].
    Δx      : Tuple or strucutre containint the horizontal and vertical grid resolution [ m ]
    Δy      : Tuple or strucutre containint the horizontal and vertical grid resolution [ m ]
    NC      : Tuple or structure containing the number of horizontal and vertical 
              centroids.
    BC      : Tuple or structure containint the boundary conditions.    
    K       : Coefficient matrix in sparse format
    rhs     : Known right-hand side vector
    Num     : Tuple or structure containing the global centroid numbering
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
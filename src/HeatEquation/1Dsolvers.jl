@doc raw"""
    SolveHeat1D_explicit!( Tnew, T_ex, κ, Δt, nc, Δx, BC )

Function to solve the 1D heat equation (diffusion only, no internal
heating, konstant thermal parameters) using an explicit (foward euler) 
finite difference scheme. The equation has the form of: 
|
| ∂T/∂t = κ ∂²T / ∂x².
| 
The temperature is defined on central nodes and the heat flux on the 
vertices. Boundary conditions are currently limited to Dirichlet and 
Neumann. Using central temperature nodes requires external ghost nodes, 
which are used to define the boundary conditions. 

Tnew    : New temperature vector [ C ]
T_ex    : Temperature vector including the ghost nodes
κ       : Diffusivity [ m²/s ]
Δt      : Time step [ s ]
nc      : Number of central nodes
Δx      : Grid spacing [ m ]
BC      : Structure for the boundary condition

"""

using ExtendableSparse

function ForwardEuler!( explicit, κ, Δt, nc, Δx, BC )
    # =================================================================== #
    # LF; 19.09.2024 - Version 1.0 - Julia                                #
    # =================================================================== #
    explicit.T_ex[2:end-1]  =   explicit.T0
    # Define boundary conditions ---------------------------------------- #
    # West ---
    explicit.T_ex[1]    =   (BC.type.W==:Dirichlet) * (2 * BC.val.W - explicit.T_ex[2]) + 
                            (BC.type.W==:Neumann) * (explicit.T_ex[2] - BC.val.W*Δx)
    # East --
    explicit.T_ex[end]  =   (BC.type.W==:Dirichlet) * (2 * BC.val.E - explicit.T_ex[nc+1]) +
                            (BC.type.W==:Neumann) * (explicit.T_ex[nc+1] + BC.val.E*Δx)
    for i = 1:nc
        # Calculate temperature at point i for the new time ---        
        explicit.T[i] =   explicit.T_ex[i+1] + κ * Δt * 
                (explicit.T_ex[i + 2] - 2.0 * explicit.T_ex[i+1] + explicit.T_ex[i]) / Δx^2
    end    
end

function BackwardEuler!( implicit, nc, Δx, κ, Δt, BC, K )
    # =================================================================== #
    # LF; 19.09.2024 - Version 1.0 - Julia                                #
    # =================================================================== #
    # Define coefficients ---
    a   =   κ / Δx^2
    b   =   1 / Δt
    # Multiply rhs with 1/Δt ---    
    implicit.T0  .=   b .* implicit.T0
    # Loop over the grid points ---
    for i = 1:nc  
        # Equation number ---
        ii          =   i
        # Stencil ---
        iW          =   ii - 1      # West
        iC          =   ii          # Central
        iE          =   ii + 1      # East
        # Boundaries ---
        # If an West index is required ---
        inW    =  i==1    ? false  : true
        DirW   = (i==1    && BC.type.W==:Dirichlet) ? 1. : 0.
        NeuW   = (i==1    && BC.type.W==:Neumann  ) ? 1. : 0.
        # If an East index is required ---
        inE    =  i==nc ? false  : true
        DirE   = (i==nc && BC.type.E==:Dirichlet) ? 1. : 0.
        NeuE   = (i==nc && BC.type.E==:Neumann  ) ? 1. : 0.
        if inE K[ii,iE]    = - a end
        K[ii,iC]        =   (2 + DirW + DirE - NeuW - NeuE)*a + b
        if inW K[ii,iW]    = - a end            
        # Modify right hand side due to boundary conditions ------------- #
        if inE 
            # West boundary ---
            implicit.T0[1]   = implicit.T0[1] + 2*a*BC.val.W * DirW - a*BC.val.W*Δx * NeuW
                    end
        if inW
            # East boundary ---
            implicit.T0[nc]  = implicit.T0[nc] + 2*a*BC.val.E * DirE + a*BC.val.E*Δx * NeuE
        end
    end
    # ------------------------------------------------------------------- #    
    # Calculate temperature at new time step ---------------------------- #
    implicit.T  .=   K \ implicit.T0
    # ------------------------------------------------------------------- #    
end

function CNV!( cnv, nc, κ, Δt, Δx, BC, K1, K2 )
# ======================================================================= #
# LF; 19.09.2024 - Version 1.0 - Julia                                    #
# ======================================================================= #    
    rhs     = zeros(length(cnv.T0))
    # Define coefficients ---
    a       =   κ / 2 / Δx^2
    b       =   1 / Δt
    # Loop over the grid points ---
    for i = 1:nc
        # Equation number ---
        ii          =   i
        # Stencil ---
        iW          =   ii - 1      # West
        iC          =   ii          # Central
        iE          =   ii + 1      # East
        # Boundaries ---
        # If an West index is required ---
        inW    =  i==1    ? false  : true
        DirW   = (i==1    && BC.type.W==:Dirichlet) ? 1. : 0.
        NeuW   = (i==1    && BC.type.W==:Neumann  ) ? 1. : 0.
        # If an East index is required ---
        inE    =  i==nc ? false  : true
        DirE   = (i==nc && BC.type.E==:Dirichlet) ? 1. : 0.
        NeuE   = (i==nc && BC.type.E==:Neumann  ) ? 1. : 0.
        if inE
            K1[ii,iE]   =   - a
            K2[ii,iE]   =   a            
        end
        K1[ii,iC]       =   b + (2 + DirW + DirE - NeuW - NeuE) * a
        K2[ii,iC]       =   b - (2 + DirW + DirE - NeuW - NeuE) * a 
        if inW 
            K1[ii,iW]   =   - a
            K2[ii,iW]   =   a
        end                    
    end
    # ------------------------------------------------------------------- #
    # Berechnung der rechten Seite -------------------------------------- #
    rhs     .=   K2 * cnv.T0 
    # ------------------------------------------------------------------- #        
    # Aenderung der rechten Seite durch die Randbedingungen ------------- #    
    for i = 1:nc        
        # Boundaries         
        DirW   = (i==1    && BC.type.W==:Dirichlet) ? 1. : 0.
        NeuW   = (i==1    && BC.type.W==:Neumann  ) ? 1. : 0.
        DirE   = (i==nc && BC.type.E==:Dirichlet) ? 1. : 0.
        NeuE   = (i==nc && BC.type.E==:Neumann  ) ? 1. : 0.
        # West boundary
        rhs[1]   = rhs[1] + 4*a*BC.val.W * DirW - 2*a*BC.val.W*Δx * NeuW        
    
        # East boundary
        rhs[nc]  = rhs[nc] + 4*a*BC.val.E * DirE + 2*a*BC.val.E*Δx * NeuE
    end
    # ------------------------------------------------------------------- #    
    # Compute new temperature ------------------------------------------- #
    cnv.T      .=    K1 \ rhs
    # ------------------------------------------------------------------- #
end

function ComputeResiduals!( dc, BC, κ, Δx, Δt )
    #ComputeResiduals!(R, T, T_ex, Told, ∂T2∂x2, BC, κ, Δx, Δt)    
    # Assign temperature to extra field --------------------------------- #
    dc.T_ex[2:end-1]    .=   dc.T    
    # Define temperature on the ghost nodes; West 
    dc.T_ex[1]          =   (BC.type.W==:Dirichlet)*(2*BC.val.W - dc.T_ex[2]) + 
                            (BC.type.W==:Neumannn)*(dc.T_ex[2] - BC.val.W*Δx)
    # Define temperature on the ghost nodes; East ---
    dc.T_ex[end]        =   (BC.type.W==:Dirichlet)*(2*BC.val.E - dc.T_ex[end-1]) + 
                            (BC.type.W==:Neumannn)*(dc.T_ex[end-1] + BC.val.E*Δx)
    # ------------------------------------------------------------------- #
    # Calculate temperature derivative ---------------------------------- #
    @. dc.∂T2∂x2        =   κ * 
            (dc.T_ex[3:end] - 2 * dc.T_ex[2:end-1] + dc.T_ex[1:end-2]) / Δx^2
    # ------------------------------------------------------------------- #
    # Calculate residual ------------------------------------------------ #
    @. dc.R             =   (dc.T - dc.T0)/Δt - dc.∂T2∂x2
    # ------------------------------------------------------------------- #
end

function AssembleMatrix!( K, BC, nc, κ, Δx, Δt )
    # Define coefficients ---
    a   =   κ / Δx^2
    b   =   1 / Δt
    # Loop over the grid points ---    
    for i = 1:nc  
        # Equation number ---
        ii          =   i
        # Stencil ---
        iW          =   ii - 1      # West
        iC          =   ii          # Central
        iE          =   ii + 1      # East
        # Boundaries ---
        # If an West index is required ---
        inW    =  i==1    ? false  : true
        DirW   = (i==1    && BC.type.W==:Dirichlet) ? 1. : 0.
        NeuW   = (i==1    && BC.type.W==:Neumann  ) ? 1. : 0.
        # If an East index is required ---
        inE    =  i==nc ? false  : true
        DirE   = (i==nc && BC.type.E==:Dirichlet) ? 1. : 0.
        NeuE   = (i==nc && BC.type.E==:Neumann  ) ? 1. : 0.
        if inE
            K[ii,iE]    = - a
        end
        K[ii,iC]        =   (2 + DirW + DirE - NeuW - NeuE)*a + b
        if inW 
            K[ii,iW]    = - a
        end        
    end
    flush!(K)
end
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

function ForwardEuler!( Tnew, T_ex, κ, Δt, nc, Δx, BC )
    # =================================================================== #
    # LF; 19.09.2024 - Version 1.0 - Julia                                #
    # =================================================================== #
    # Define boundary conditions ---------------------------------------- #
    # West ---
    T_ex[1]     =   (BC.type.W==:Dirichlet) * (2 * BC.val.W - T_ex[2]) + 
                    (BC.type.W==:Neumann) * (T_ex[2] - BC.val.W*Δx)
    # East --
    T_ex[end]   =   (BC.type.W==:Dirichlet) * (2 * BC.val.E - T_ex[nc+1]) +
                    (BC.type.W==:Neumann) * (T_ex[nc+1] + BC.val.E*Δx)
    for i = 1:nc
        # Calculate temperature at point i for the new time ---        
        Tnew[i] =   T_ex[i+1] + κ * Δt * 
                            (T_ex[i + 2] - 2.0 * T_ex[i+1] + T_ex[i]) / Δx^2
    end
    return Tnew
end

function BackwardEuler!( T0, nc, Δx, κ, Δt, BC, K )
    # =================================================================== #
    # LF; 19.09.2024 - Version 1.0 - Julia                                #
    # =================================================================== #
    # Define coefficients ---
    a   =   κ / Δx^2
    b   =   1 / Δt
    # Multiply rhs with 1/Δt ---
    T0  =   b * T0
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
            T0[1]   = T0[1] + 2*a*BC.val.W * DirW - a*BC.val.W*Δx * NeuW
                    end
        if inW
            # East boundary ---
            T0[nc]  = T0[nc] + 2*a*BC.val.E * DirE + a*BC.val.E*Δx * NeuE
        end
    end
    # ------------------------------------------------------------------- #    
    # Calculate temperature at new time step ---------------------------- #
    T1      =   K\T0
    # ------------------------------------------------------------------- #    
    return T1, K 
end

function CNV!(T0,nc,κ,Δt,Δx,BC,K1,K2)
    # ======================================================================= #
    # LF; 19.09.2024 - Version 1.0 - Julia                                    #
    # ======================================================================= #    
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
        rhs     =   K2*T0 
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
        T1      = K1\rhs
        # ------------------------------------------------------------------- #
        return T1, K1, K2
    end

function ComputeResiduals!(R, T, T_ex, Told, ∂T2∂x2, BC, κ, Δx, Δt)    
    # Assign temperature to extra field --------------------------------- #
    T_ex[2:end-1]       =   T    
    # Define temperature on the ghost nodes; West 
    T_ex[1]             =   (BC.type.W==:Dirichlet)*(2*BC.val.W - T_ex[2]) + 
                            (BC.type.W==:Neumannn)*(T_ex[2] - BC.val.W*Δx)
    # Define temperature on the ghost nodes; East ---
    T_ex[end]           =   (BC.type.W==:Dirichlet)*(2*BC.val.E - T_ex[end-1]) + 
                            (BC.type.W==:Neumannn)*(T_ex[end-1] + BC.val.E*Δx)
    # ------------------------------------------------------------------- #
    # Calculate temperature derivative ---------------------------------- #
    ∂T2∂x2              =   κ .* 
            (T_ex[3:end] .- 2 .* T_ex[2:end-1] .+ T_ex[1:end-2]) ./ Δx^2
    # ------------------------------------------------------------------- #
    # Calculate residual ------------------------------------------------ #
    R                   =   (T .- Told)./Δt .- ∂T2∂x2
    # ------------------------------------------------------------------- #
end

function AssembleMatrix!(K,BC,nc,κ,Δx,Δt)
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

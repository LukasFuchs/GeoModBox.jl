using ExtendableSparse
@doc raw"""
    ForwardEuler1Dc!( explicit, κ, Δx, Δt, nc, BC)

Solves the onedimensional heat diffusion equation assuming no internal heating and
constant thermal parameters using an explicit, forward euler finite difference scheme.

The temperature is defined on central nodes and the heat flux on the vertices. 
Boundary conditions are currently limited to Dirichlet and Neumann. Using central 
temperature nodes requires external ghost nodes, which are used to define the 
boundary conditions. 

    explicit    : Tuple, containing the regular temperature array T and 
                  array containing the ghost nodes T_ex
    κ           : Diffusivity [ m²/s ]
    Δt          : Time step [ s ]
    nc          : Number of central nodes
    Δx          : Grid spacing [ m ]
    BC          : Tuple for the boundary condition
"""
function ForwardEuler1Dc!( explicit, κ, Δx, Δt, nc, BC )
    # =================================================================== #
    # LF; 19.09.2024 - Version 1.0 - Julia                                #
    # =================================================================== #
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
    explicit.T_ex[2:end-1]  .=  explicit.T
end

@doc raw"""
    ForwardEuler1D!( explicit, κ, Δx, Δt, nc, BC)

Solves the onedimensional heat diffusion equation assuming internal heating and
variable thermal parameters using an explicit, forward euler finite difference scheme.

The temperature is defined on central nodes and the heat flux on the vertices. 
Boundary conditions are currently limited to Dirichlet and Neumann. Using central 
temperature nodes requires external ghost nodes, which are used to define the 
boundary conditions. 

    T           : Tuple, containing the regular temperature array T and 
                  array containing the ghost nodes T_ex
    Py          : Tuple, containing the thermal parameters ρ, k, cp, and H [ W/kg ]
    Δt          : Time step [ s ]
    Δy          : Grid spacing [ m ]
    nc          : Number of central nodes
    BC          : Tuple for the boundary condition

"""
function ForwardEuler1D!(T,Py,Δt,Δy,nc,BC)
    
    if size(Py.k,1) == 1
        k   =   Py.k.*ones(nc+1,1)
        ρ   =   Py.ρ.*ones(nc,1)
        cp  =   Py.cp.*ones(nc,1)
    else
        k   =   Py.k
        ρ   =   Py.ρ
        cp  =   Py.cp
    end
    if size(Py.H,1) == 1
        H   =   Py.H.*ones(nc,1)      #   [H] = W/kg; [Q] = [ρ*H], [Q] = W/m³
    else
        H   =   Py.H
    end

    # Define boundary conditions ---------------------------------------- #
    # South ---
    T.T_ex[1]   =   (BC.type.S==:Dirichlet) * (2 * BC.val.S - T.T_ex[2]) + 
                    (BC.type.S==:Neumann) * (T.T_ex[2] - BC.val.S*Δy)
    # North ---
    T.T_ex[end] =   (BC.type.N==:Dirichlet) * (2 * BC.val.N - T.T_ex[nc+1]) +
                    (BC.type.N==:Neumann) * (T.T_ex[nc+1] + BC.val.N*Δy)
    
    for j = 1:nc
        a       =   Δt/(Δy^2*ρ[j]*cp[j])
        T.T[j]  =   a*k[j]*T.T_ex[j] + 
                    (1-a*(k[j+1]+k[j]))*T.T_ex[j+1] +
                    a*k[j+1]*T.T_ex[j+2] +
                    H[j]*Δt/cp[j]
    end
    T.T_ex[2:end-1]     .= T.T
end

@doc raw"""
    BackwardEuler1Dc!( implicit, κ, Δx, Δt, nc, BC , K)

Solves the onedimensional heat diffusion equation assuming no internal heating and
constant thermal parameters using an implicit, backward euler finite difference scheme.

The temperature is defined on central nodes and the heat flux on the vertices. 
Boundary conditions are currently limited to Dirichlet and Neumann. Using central 
temperature nodes requires external ghost nodes, which are used to define the 
boundary conditions. 

    implicit    : Tuple, containing the current temperature array T0 and 
                  the new temperature array T
    κ           : Diffusivity [ m²/s ]
    Δt          : Time step [ s ]
    nc          : Number of central nodes
    Δx          : Grid spacing [ m ]
    BC          : Tuple for the boundary condition
    K           : Coefficient matrix for linear system of equations
"""
function BackwardEuler1Dc!( implicit, κ, Δx, Δt, nc, BC , K, rhs)
    # =================================================================== #
    # LF; 19.09.2024 - Version 1.0 - Julia                                #
    # =================================================================== #
    # Define coefficients ---
    a   =   κ / Δx^2
    b   =   1 / Δt
    # Multiply rhs with 1/Δt ---    
    rhs  .=   b .* implicit.T
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
        # West boundary ---
        rhs[i]  +=  2*a*BC.val.W * DirW - 
                        a*BC.val.W*Δx * NeuW +
                        2*a*BC.val.E * DirE + 
                        a*BC.val.E*Δx * NeuE
    end
    # ------------------------------------------------------------------- #    
    # Calculate temperature at new time step ---------------------------- #
    implicit.T  .=   K \ rhs
    # ------------------------------------------------------------------- #    
end

@doc raw"""
    CNA1Dc!( cna, κ, Δx, Δt, nc, BC, K1, K2 )

Solves the onedimensional heat diffusion equation assuming no internal heating and
constant thermal parameters using Crank-Nicolson finite difference scheme.

The temperature is defined on central nodes and the heat flux on the vertices. 
Boundary conditions are currently limited to Dirichlet and Neumann. Using central 
temperature nodes requires external ghost nodes, which are used to define the 
boundary conditions. 

    cna         : Tuple, containing the current temperature array T0 and 
                  the new temperature array T
    κ           : Diffusivity [ m²/s ]
    Δt          : Time step [ s ]
    nc          : Number of central nodes
    Δx          : Grid spacing [ m ]
    BC          : Tuple for the boundary condition
    K1          : Coefficient matrix for the unknow variables 
    K2          : Coefficient matrix for the know variables
"""
function CNA1Dc!( cna, κ, Δx, Δt, nc, BC, K1, K2 )
# ======================================================================= #
# LF; 19.09.2024 - Version 1.0 - Julia                                    #
# ======================================================================= #    
    rhs     = zeros(length(cna.T))
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
    rhs     .=   K2 * cna.T
    # ------------------------------------------------------------------- #        
    # Aenderung der rechten Seite durch die Randbedingungen ------------- #    
    for i = 1:nc        
        # Boundaries         
        DirW   = (i==1    && BC.type.W==:Dirichlet) ? 1. : 0.
        NeuW   = (i==1    && BC.type.W==:Neumann  ) ? 1. : 0.
        DirE   = (i==nc && BC.type.E==:Dirichlet) ? 1. : 0.
        NeuE   = (i==nc && BC.type.E==:Neumann  ) ? 1. : 0.
        # West boundary
        rhs[i]  +=  4*a*BC.val.W * DirW - 
                        2*a*BC.val.W*Δx * NeuW +
                        4*a*BC.val.E * DirE + 
                        2*a*BC.val.E*Δx * NeuE
    end
    # ------------------------------------------------------------------- #    
    # Compute new temperature ------------------------------------------- #
    cna.T      .=    K1 \ rhs
    # ------------------------------------------------------------------- #
end

@doc raw"""
    ComputeResiduals1Dc!( cna, κ, Δx, Δt, nc, BC, K1, K2 )

Computes the residual of the onedimensional heat diffusion equation assuming 
no internal heating and constant thermal parameters.

The temperature is defined on central nodes and the heat flux on the vertices. 
Boundary conditions are currently limited to Dirichlet and Neumann. Using central 
temperature nodes requires external ghost nodes, which are used to define the 
boundary conditions. 

    dc          : Tuple, containing the current temperature array T, 
                  the temperature array with ghost nodes T_ex,
                  the partial derivatives ∂T2∂x2, and the
                  residual R
    κ           : Diffusivity [ m²/s ]
    Δx          : Grid spacing [ m ]
    Δt          : Time step [ s ]       
    BC          : Tuple for the boundary condition
"""
function ComputeResiduals1Dc!( dc, κ, Δx, Δt, BC )
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

@doc raw"""
    AssembleMatrix1Dc!( κ, Δx, Δt, nc, BC, K )

Setup the coefficient matrix for the linear system of equations. 
    
"""
function AssembleMatrix1Dc!( κ, Δx, Δt, nc, BC, K )
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
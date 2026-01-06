using ExtendableSparse
# ======================================================================= #
# Time-dependent solvers, constant thermal parameters =================== #
# ======================================================================= #
"""
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

"""
    ComputeResiduals1Dc!( cna, κ, Δx, Δt, nc, BC, K1, K2 )

Computes the residual of the onedimensional heat diffusion equation assuming 
no internal heating and constant thermal parameters.

The temperature is defined on central nodes and the heat flux on the vertices. 
Boundary conditions are currently limited to Dirichlet and Neumann. Using central 
temperature nodes requires external ghost nodes, which are used to define the 
boundary conditions. 

    dc          : Tuple, containing the current temperature array T, 
                  the temperature array with ghost nodes T_ex,
                  the partial derivatives ∂2T∂x2, and the
                  residual R
    κ           : Diffusivity [ m²/s ]
    Δx          : Grid spacing [ m ]
    Δt          : Time step [ s ]       
    BC          : Tuple for the boundary condition
"""
function ComputeResiduals1Dc!( R, T, T_ex, T0, T_ex0, ∂2T, Q, ρ, cp, κ, BC, Δx, Δt;C=0)
    if C < 1
        # Implicit
        T_ex[2:end-1]   .=  T    
        T_ex[1]         =   (BC.type.W==:Dirichlet)*(2*BC.val.W - T_ex[2]) + (BC.type.W==:Neumannn)*(T_ex[2] - BC.val.W*Δx)
        T_ex[end]       =   (BC.type.W==:Dirichlet)*(2*BC.val.E - T_ex[end-1]) + (BC.type.W==:Neumannn)*(T_ex[end-1] + BC.val.E*Δx)
        @. ∂2T.∂x2      =   (T_ex[3:end] - 2 * T_ex[2:end-1] + T_ex[1:end-2]) / Δx^2
        if C==0.5
            # CNA
            T_ex0[2:end-1]  .=  T0    
            T_ex0[1]        =   (BC.type.W==:Dirichlet)*(2*BC.val.W - T_ex0[2]) + (BC.type.W==:Neumannn)*(T_ex0[2] - BC.val.W*Δx)
            T_ex0[end]      =   (BC.type.W==:Dirichlet)*(2*BC.val.E - T_ex0[end-1]) + (BC.type.W==:Neumannn)*(T_ex0[end-1] + BC.val.E*Δx)
            @. ∂2T.∂x20     =   (T_ex0[3:end] - 2 * T_ex0[2:end-1] + T_ex0[1:end-2]) / Δx^2
        end
    else
        # explicit
        T_ex0[2:end-1]  .=   T0    
        T_ex0[1]        =   (BC.type.W==:Dirichlet)*(2*BC.val.W - T_ex0[2]) + 
                                (BC.type.W==:Neumannn)*(T_ex0[2] - BC.val.W*Δx)
        T_ex0[end]      =   (BC.type.W==:Dirichlet)*(2*BC.val.E - T_ex0[end-1]) + 
                                (BC.type.W==:Neumannn)*(T_ex0[end-1] + BC.val.E*Δx)
        @. ∂2T.∂x20     =   (T_ex0[3:end] - 2 * T_ex0[2:end-1] + T_ex0[1:end-2]) / Δx^2
    end
    # Calculate residual ------------------------------------------------ #
    @. R             =   (T - T0)/Δt - κ*((1-C)*∂2T.∂x2 + C*∂2T.∂x20) - Q/ρ/cp
    # ------------------------------------------------------------------- #
end

"""
    AssembleMatrix1Dc!( κ, Δx, Δt, nc, BC, K; C )

Setup the coefficient matrix for the linear system of equations. 
    
"""
function AssembleMatrix1Dc!( κ, Δx, Δt, nc, BC, K;C=0 )
    # Define coefficients ---
    a   =   (κ*(1-C)) / Δx^2
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



"""
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

"""
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
# ======================================================================= #
# Time-dependent solvers, variable thermal parameters =================== #
# ======================================================================= #
"""
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

"""
    ComputeResiduals1D!(R, T, T_ex, T0, T_ex0, Q, ∂T, q, ρ, Cp, k, BC, Δx, Δt;C=0)

Function to calculate the residual of the one dimensional heat diffusion equation assuming 
variable thermal parameters and radiogenic heating only. The residual is required for 
defection correction to solve the linear system of equations. 

The temperature is defined on central nodes and the heat flux on the vertices. 
Boundary conditions are currently limited to Dirichlet and Neumann. Using central 
temperature nodes requires external ghost nodes, which are used to define the 
boundary conditions. The caluclated residual is used in an iteration loop to calculate 
the correction term of the initial temperature guess. The coefficient matrix is build 
via `AssembleMatrix1D()`. 

    R           : 1D array for the residual
    T           : 1D centroid temperature field of the next time step
    T_ex        : 1D extended centroid temperature field including the ghost nodes
    T_0         : 1D centroid temperature field of the current time step
    T_ex0       : 1D extended centroid temperature field of the current time step
    Q           : Volumetric heat production rate [ W/m³ ]
    ∂T          : Tuple containing the first space derivatives of the temperature
                  in the horizontal direction
    q           : Tuple or structure containing the 2D horizontal heat flux
    ρ           : Density [ kg/m³ ]
    Cp          : Specific heat capacity [ J/kg/K ]
    k           : Thermal conductivity [ W/m/K ]
    BC          : Tuple for the boundary condition
    Δx          : Tuple or strucutre containint the horizontal and vertical grid resolution [ m ]
    Δt          : Time step [ s ]

Optional input values: 
    C           : Constant defining the residual for a certain finite difference discretization 
                  method, i.e.: 
                        C = 0   -> implicit, backward Euler discretization (default)
                        C = 0.5 -> Crank-Nicolson discretization
                        C = 1   -> explicit, forward Euler discretization
    
"""
function ComputeResiduals1D!(R, T, T_ex, T0, T_ex0, Q, ∂T, q, ρ, Cp, k, BC, Δx, Δt;C=0)
    if C < 1
        @. T_ex[2:end-1]    = T 
        T_ex[1]          = (BC.type.W==:Dirichlet) * (2*BC.val.W - T_ex[    2]) + (BC.type.W==:Neumann) * (T_ex[    2] - Δx/k[  1]*BC.val.W)
        T_ex[end]        = (BC.type.E==:Dirichlet) * (2*BC.val.E - T_ex[end-1]) + (BC.type.E==:Neumann) * (T_ex[end-1] + Δx/k[end]*BC.val.E)
        @. ∂T.∂x = (T_ex[2:end] - T_ex[1:end-1])/Δx
        @. q.x   = -k * ∂T.∂x
        if C==0.5
            @. T_ex0[2:end-1]   = T0 
            T_ex0[1]         = (BC.type.W==:Dirichlet) * (2*BC.val.W - T_ex0[    2]) + (BC.type.W==:Neumann) * (T_ex0[    2] - Δx/k[  1]*BC.val.W)
            T_ex0[end]       = (BC.type.E==:Dirichlet) * (2*BC.val.E - T_ex0[end-1]) + (BC.type.E==:Neumann) * (T_ex0[end-1] + Δx/k[end]*BC.val.E)
            @. ∂T.∂x0   =   (T_ex0[2:end] - T_ex0[1:end-1])/Δx
            @. q.x0     =   -k * ∂T.∂x0
        end
    else
        @. T_ex0[2:end-1]   = T0 
        T_ex0[1]         = (BC.type.W==:Dirichlet) * (2*BC.val.W - T_ex0[    2]) + (BC.type.W==:Neumann) * (T_ex0[    2] - Δx/k[  1]*BC.val.W)
        T_ex0[end]       = (BC.type.E==:Dirichlet) * (2*BC.val.E - T_ex0[end-1]) + (BC.type.E==:Neumann) * (T_ex0[end-1] + Δx/k[end]*BC.val.E)
        @. ∂T.∂x0   =   (T_ex0[2:end] - T_ex0[1:end-1])/Δx
        @. q.x0     =   -k * ∂T.∂x0
        @. q.x      =   q.x0
    end
    @. R     = ρ*Cp*(T - T0)/Δt + 
                    (1-C)*((q.x[2:end] - q.x[1:end-1])/Δx) + 
                    C*((q.x0[2:end] - q.x0[1:end-1])/Δx) - 
                    Q
end

"""
    AssembleMatrix1D(ρ, cp, k, Δx, Δt, nc, BC, K;C=0 )

Function to build the coefficient matrix K for the unknown centroid temperature field in the 
2D heat diffusion equation assuming variable thermal parameter and radiogenig heating only. 

The coefficient matrix is build using a five-point finite difference stencil, resulting in 
a five, non-zero diagonal matrix to solve the system of equations. 

    ρ       : Density [ kg/m³ ]
    cp      : Specific heat capacity [ J/kg/K ]
    k       : Thermal conductivity [ W/m/K ]
    Δx      : Tuple or structure containing the horizontal grid resolution
    Δt      : Time step [ s ]
    nc      : Tuple or structure containing the number of centroids in the horizontal direction 
    BC      : Tuple for the boundary condition
    K       : Coefficient matrix

Optional input values: 
    C           : Constant defining the residual for a certain finite difference discretization 
                  method, i.e.: 
                        C = 0   -> implicit, backward Euler discretization (default)
                        C = 0.5 -> Crank-Nicolson discretization
                        C = 1   -> explicit, forward Euler discretization
"""
function AssembleMatrix1D( ρ, cp, k, Δx, Δt, nc, BC, K;C=0 )
    for i=1:nc
        # Equation number
        ii = i
        # Stencil
        iW = ii - 1
        iC = ii
        iE = ii + 1
        # Boundaries
        inW    =  i==1    ? false  : true   
        DirW   = (i==1    && BC.type.W==:Dirichlet) ? 1. : 0.
        NeuW   = (i==1    && BC.type.W==:Neumann  ) ? 1. : 0.
        inE    =  i==nc ? false  : true   
        DirE   = (i==nc && BC.type.E==:Dirichlet) ? 1. : 0.
        NeuE   = (i==nc && BC.type.E==:Neumann  ) ? 1. : 0.
        # Material coefficient
        kC = k[i]*(1-C)
        kE = k[i+1]*(1-C)
        # Linear system coefficients
        if inW K[ii,iW] = kC .* (DirW + NeuW - 1) ./ Δx .^ 2 end
        K[ii,iC] = cp[i] .* ρ[i] ./ Δt + (-kE .* (-DirE + NeuE - 1) ./ Δx + kC .* (DirW - NeuW + 1) ./ Δx) ./ Δx
        if inE K[ii,iE] = -kE .* (-DirE - NeuE + 1) ./ Δx .^ 2 end
    end
    return flush!(K)
end


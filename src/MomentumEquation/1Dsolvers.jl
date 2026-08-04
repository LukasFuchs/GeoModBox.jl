using ExtendableSparse

"""
    ComputeStokesResiduals1D!(D, ∂P∂x, Δy, BC)

Compute the residual of the one-dimensional incompressible Stokes equation

    -∂P∂x + ∂τxy/∂y = 0,

where the shear stress is given by

    τxy = η ∂vₓ/∂y.

The function applies the prescribed Dirichlet or Neumann boundary
conditions using ghost nodes, computes the velocity gradient,
shear stress, and stress divergence, and stores the resulting
residual in `D.R`.

# Arguments
    - `D`   : Structure containing the velocity field, viscosity, stress,
                derivatives, and residual arrays.
    - `∂P∂x`: Prescribed pressure gradient (scalar or array).
    - `Δy`  : Grid spacing in the y-direction.
    - `BC`  : Boundary condition structure specifying the boundary types
    and values.

# Returns
Nothing. The residual and intermediate quantities are updated in-place.
"""
function ComputeStokesResiduals1D!( D, ∂P∂x, Δy, BC)
    @. D.vₓₑ[2:end-1]     =   D.vₓ
    # Bottom ---
    D.vₓₑ[1]            =   (BC.type.S==:Dirichlet)*(2*BC.val.S - D.vₓₑ[2]) + 
                            (BC.type.S==:Neumann)*(D.vₓₑ[2] - Δy*BC.val.S)  
    # Top --- 
    D.vₓₑ[end]          =   (BC.type.N==:Dirichlet)*(2*BC.val.N - D.vₓₑ[end-1]) + 
                            (BC.type.N==:Neumann)*(D.vₓₑ[end-1] + Δy*BC.val.N)  
    # ---
    @. D.∂vₓ∂y          =   (D.vₓₑ[2:end] - D.vₓₑ[1:end-1]) / Δy
    @. D.τxy            =   D.η * D.∂vₓ∂y
    @. D.∂τxy∂y         =   (D.τxy[2:end] - D.τxy[1:end-1]) / Δy
    @. D.R              =   -∂P∂x + D.∂τxy∂y
end

"""
    AssembleStokesMatrix1D(nc, η, Δy, BC, K)

Assemble the sparse coefficient matrix for the one-dimensional
incompressible Stokes equation

    ∂/∂y(η ∂vₓ/∂y) = ∂P/∂x.

The function constructs the variable-viscosity finite-difference
stencil and incorporates the prescribed Dirichlet or Neumann boundary
conditions into the matrix coefficients. A system with Neumann
conditions at both boundaries is rejected because the resulting matrix
is singular.

# Arguments
    - `nc`  : Number of computational cells.
    - `η`   : Dynamic viscosity defined at the cell faces.
    - `Δy`  : Grid spacing in the y-direction.
    - `BC`  : Boundary condition structure specifying the southern and
                northern boundary types and values.
    - `K`   : Extendable sparse matrix used to store the coefficient matrix.

# Returns
The finalized sparse coefficient matrix.
"""
function AssembleStokesMatrix1D(nc, η, Δy, BC, K)
    if BC.type.S == :Neumann && BC.type.N == :Neumann
        throw(ArgumentError(
            "The 1-D Stokes system is singular when Neumann conditions " *
            "are applied at both boundaries. Prescribe at least one " *
            "Dirichlet condition or impose an additional velocity constraint."
        ))
    end
    # Zusammenstellen der Koeffizientenmatrix ------------------------------- #
    for i = 1:nc
        a   =   η[i] / Δy^2.0
        b   =   -(η[i]+η[i+1]) / Δy^2.0
        c   =   η[i+1] / Δy^2.0
        # Gleichungsnummer ---
        ii  =   i
        # Stempel ---
        iN  =   ii + 1      #   Norden
        iC  =   ii          #   Zentral
        iS  =   ii - 1      #   Süden
        # Ränder ---
        # Falls ein Süd Index gebrauch wird ---
        inS    =  i==1    ? false  : true
        DirS   = (i==1    && BC.type.S==:Dirichlet) ? 1. : 0.
        NeuS   = (i==1    && BC.type.S==:Neumann  ) ? 1. : 0.
        # If an East index is required ---
        inN    =  i==nc ? false  : true
        DirN   = (i==nc && BC.type.N==:Dirichlet) ? 1. : 0.
        NeuN   = (i==nc && BC.type.N==:Neumann  ) ? 1. : 0.
        if inS K[ii,iS]    = a end
            K[ii,iC]       = b + (NeuS - DirS)*a + (NeuN - DirN)*c
        if inN K[ii,iN]    = c end    
    end
    return flush!(K)
end

"""
    Stokes_1D_direct(vₓ, η, Δy, nc, BC, K, rhs)

Solve the one-dimensional incompressible Stokes equation using a
direct sparse linear solver.

The function assembles the finite-difference coefficient matrix,
incorporates the prescribed Dirichlet or Neumann boundary conditions,
modifies the right-hand side accordingly, and solves the resulting
linear system for the horizontal velocity. A system with Neumann
conditions at both boundaries is rejected because the resulting matrix
is singular.

# Arguments
    - `vₓ`  : Array storing the computed horizontal velocity.
    - `η`   : Dynamic viscosity defined at the cell faces.
    - `Δy`  : Grid spacing in the y-direction.
    - `nc`  : Number of computational cells.
    - `BC`  : Boundary condition structure specifying the southern and
                northern boundary types and values.
    - `K`   : Extendable sparse matrix used to assemble the coefficient matrix.
    - `rhs` : Right-hand side vector representing the pressure gradient and
    any additional body forces.

# Returns
The computed horizontal velocity `vₓ`.
"""
function Stokes_1D_direct(vₓ,η,Δy,nc,BC,K,rhs)

    if BC.type.S == :Neumann && BC.type.N == :Neumann
        throw(ArgumentError(
            "The 1-D Stokes system is singular when Neumann conditions " *
            "are applied at both boundaries. Prescribe at least one " *
            "Dirichlet condition or impose an additional velocity constraint."
        ))
    end

    rhs_bc = copy(rhs)

    # Zusammenstellen der Koeffizientenmatrix ------------------------------- #
    for i = 1:nc
        a   =   η[i] / Δy^2.0
        b   =   -(η[i]+η[i+1]) / Δy^2.0
        c   =   η[i+1] / Δy^2.0
        # Gleichungsnummer ---
        ii  =   i
        # Stempel ---
        iN  =   ii + 1      #   Norden
        iC  =   ii          #   Zentral
        iS  =   ii - 1      #   Süden
        # Ränder ---
        # Falls ein Süd Index gebrauch wird ---
        inS    =  i==1    ? false  : true
        DirS   = (i==1    && BC.type.S==:Dirichlet) ? 1. : 0.
        NeuS   = (i==1    && BC.type.S==:Neumann  ) ? 1. : 0.
        # If an East index is required ---
        inN    =  i==nc ? false  : true
        DirN   = (i==nc && BC.type.N==:Dirichlet) ? 1. : 0.
        NeuN   = (i==nc && BC.type.N==:Neumann  ) ? 1. : 0.
        if inS K[ii,iS]    = a end
            K[ii,iC]       = b + (NeuS - DirS)*a + (NeuN - DirN)*c
        if inN K[ii,iN]    = c end    
        # Änderung der rechten Seite ---
        # rhs[i]      +=  a*BC.val.S*Δy * NeuS - 
        #                     2*a*BC.val.S * DirS - 
        #                     c*BC.val.N*Δy * NeuN - 
        #                     2*c*BC.val.N * DirN
        rhs_bc[i]   +=  a * BC.val.S * Δy * NeuS -
                            2 * a * BC.val.S * DirS -
                            c * BC.val.N * Δy * NeuN -
                            2 * c * BC.val.N * DirN
    end
    # ----------------------------------------------------------------------- #
    # Lösung des Gleichungssystems ------------------------------------------ #
    # vₓ      .=   K \ rhs
    vₓ      .=  K \ rhs_bc
    # ----------------------------------------------------------------------- #
    return vₓ

end
using ExtendableSparse

# ======================================================================= #
# Time-dependent solvers, constant thermal parameters =================== #
# ======================================================================= #
"""
    ForwardEuler2Dc!(
        D, κ, Δx, Δy, Δt, NC, BC;
        Q=zeros(NC...), ρ=3300.0, cp=1200.0, Qₛ=zeros(NC...)
    )

Solves the two-dimensional transient heat equation using an explicit Forward
Euler finite-difference scheme with constant thermal diffusivity.

Temperature is defined at the cell centroids, while the diffusion term is
discretized using second-order central finite differences. Ghost nodes are
used to impose Dirichlet and Neumann boundary conditions at the western,
eastern, southern, and northern boundaries.

Optional volumetric and shear-heating source terms may be included.

# Arguments

    D       : Structure or tuple containing:
              `T`, the temperature field defined at the cell centroids, and
              `T_ex`, the extended temperature field including ghost nodes.
              Both arrays are updated in place.
    κ       : Thermal diffusivity.
    Δx      : Horizontal grid spacing.
    Δy      : Vertical grid spacing.
    Δt      : Time step.
    NC      : Structure or tuple containing the number of centroid nodes
              in the x- and y-directions.
    BC      : Structure or tuple defining the boundary-condition types and
              values at the western, eastern, southern, and northern
              boundaries.

# Keyword Arguments

    Q       : Volumetric heat-production rate (default: `zeros(NC...)`).
    ρ       : Density (default: `3300.0`).
    cp      : Specific heat capacity (default: `1200.0`).
    Qₛ      : Volumetric shear-heating rate (default: `zeros(NC...)`).

# Notes

The method is first-order accurate in time and second-order accurate in
space. Since the temperature update is explicit, the timestep must satisfy
the diffusive stability condition in two dimensions.

For a uniform grid (`Δx = Δy = h`), the stability criterion is

    Δt ≤ h² / (4κ).

For non-uniform grids, the timestep should satisfy

    Δt ≤ 1 / (2κ(1/Δx² + 1/Δy²)).

After each timestep, the updated centroid temperatures are copied to the
interior of `D.T_ex`. The ghost-node values are updated again when the
routine is called for the next timestep.

"""
function ForwardEuler2Dc!(D, κ, Δx, Δy, Δt, NC, BC; 
                Q = zeros(NC...), ρ = 3300.0, cp = 1200.0, Qₛ = zeros(NC...) )
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
                (Q[i,j] + Qₛ[i,j]) * Δt / ρ / cp 
        end
    end
    # ------------------------------------------------------------------- #
    # Update extended temperature --------------------------------------- #
    D.T_ex[2:end-1,2:end-1]     .=  D.T
    # ------------------------------------------------------------------- #    
end

"""
    ComputeResiduals2Dc!(
        R, T, T_ex, T0, T_ex0, ∂2T, κ, BC, Δ, Δt;
        C=0.0, Q=0.0, ρ=3300.0, cp=1200.0, Qₛ=0.0
    )

Computes the residual of the two-dimensional transient heat equation with
constant thermal diffusivity and optional volumetric heat sources.

Temperature is defined at the cell centroids, and the diffusion operator is
discretized using second-order central finite differences. Ghost nodes are
used to impose Dirichlet and Neumann boundary conditions at the western,
eastern, southern, and northern boundaries.

The temporal discretization is controlled by `C` and can represent Backward
Euler, Crank–Nicolson, or Forward Euler time integration. The residual is
intended for defect-correction iterations and is used together with the
coefficient matrix assembled by `AssembleMatrix2Dc`.

# Arguments

    R       : Residual field defined at the cell centroids.
    T       : Temperature field at the current iteration or new time level.
    T_ex    : Extended current-temperature field including ghost nodes.
    T0      : Temperature field at the previous time level.
    T_ex0   : Extended previous-temperature field including ghost nodes.
    ∂2T     : Structure or tuple containing the second spatial derivatives
              `∂x2`, `∂y2`, `∂x20`, and `∂y20`.
    κ       : Thermal diffusivity.
    BC      : Structure or tuple defining the boundary-condition types and
              values at the western, eastern, southern, and northern
              boundaries.
    Δ       : Structure or tuple containing the horizontal and vertical grid
              spacings `x` and `y`.
    Δt      : Time step.

# Keyword Arguments

    C       : Temporal weighting parameter (default: `0.0`):
              `0.0` for Backward Euler,
              `0.5` for Crank–Nicolson,
              `1.0` for Forward Euler.
    Q       : Volumetric heat-production rate (default: `0.0`).
    ρ       : Density (default: `3300.0`).
    cp      : Specific heat capacity (default: `1200.0`).
    Qₛ      : Volumetric shear-heating rate (default: `0.0`).

# Notes

The residual is evaluated as the sum of

- the transient temperature-change term,
- the weighted two-dimensional diffusion term, and
- the volumetric heat-source terms.

Backward Euler evaluates the diffusion operator entirely at the new time
level, Crank–Nicolson averages the previous and current time levels, and
Forward Euler uses only the previous time level.

The routine updates the ghost-node values and second spatial derivatives
before evaluating the residual.
    
"""
function ComputeResiduals2Dc!( R, T, T_ex, T0, T_ex0, ∂2T, κ, BC, Δ, Δt;
                C = 0, Q = 0.0, ρ = 3300.0, cp = 1200.0, Qₛ = 0.0 )
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
    @. R     = (T - T0)/Δt - κ*((1-C)*(∂2T.∂x2 + ∂2T.∂y2) + C*(∂2T.∂x20 + ∂2T.∂y20)) - Q/ρ/cp - Qₛ/ρ/cp
end

"""
    AssembleMatrix2Dc(κ, BC, Num, nc, Δ, Δt; C=0.0)

Assembles the coefficient matrix for the two-dimensional transient heat
equation with constant thermal diffusivity.

Temperature is defined at the cell centroids, and the diffusion operator is
discretized using second-order central finite differences. Dirichlet and
Neumann boundary conditions are incorporated directly into the matrix
coefficients at the western, eastern, southern, and northern boundaries.

The temporal discretization is controlled by `C`. For `C < 1`, the matrix
contains the implicit contribution of the diffusion operator. For `C = 1`,
the diffusion contribution vanishes and the matrix contains only the
transient storage term.

# Arguments

    κ       : Thermal diffusivity.
    BC      : Structure or tuple defining the boundary-condition types and
              values at the western, eastern, southern, and northern
              boundaries.
    Num     : Structure or tuple containing the global numbering of the
              centroid nodes.
    nc      : Structure or tuple containing the number of centroid nodes
              in the x- and y-directions.
    Δ       : Structure or tuple containing the horizontal and vertical grid
              spacings `x` and `y`.
    Δt      : Time step.

# Keyword Arguments

    C       : Temporal weighting parameter (default: `0.0`):
              `0.0` for Backward Euler,
              `0.5` for Crank–Nicolson,
              `1.0` for Forward Euler.

# Notes

The matrix corresponds to the current-time contribution of the generalized
time-discretization scheme used by `ComputeResiduals2Dc!`.

The diffusion operator is discretized using the classical five-point finite
difference stencil, resulting in a sparse matrix with five non-zero diagonals.
Boundary conditions modify the coefficients of the equations adjacent to the
domain boundaries.

The matrix is assembled internally as an `ExtendableSparseMatrix` and
finalized using `flush!(K)` before being returned.
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
    BackwardEuler2Dc!(
        D, κ, Δx, Δy, Δt, NC, BC, rhs, K, Num;
        Q=zeros(NC...), ρ=3300.0, cp=1200.0, Qₛ=zeros(NC...)
    )

Solves the two-dimensional transient heat equation using an implicit Backward
Euler finite-difference scheme with constant thermal diffusivity.

Temperature is defined at the cell centroids, and the diffusion operator is
discretized using second-order central finite differences. Dirichlet and
Neumann boundary conditions are incorporated directly into the coefficient
matrix and right-hand-side vector at the western, eastern, southern, and
northern boundaries.

Optional volumetric heat-production and shear-heating terms may be included.
The resulting linear system is solved using a direct left-matrix division.

# Arguments

    D       : Structure or tuple containing:
              `T`, the temperature field defined at the cell centroids, and
              `T_ex`, the extended temperature field including ghost nodes.
              Both arrays are updated in place.
    κ       : Thermal diffusivity.
    Δx      : Horizontal grid spacing.
    Δy      : Vertical grid spacing.
    Δt      : Time step.
    NC      : Structure or tuple containing the number of centroid nodes
              in the x- and y-directions.
    BC      : Structure or tuple defining the boundary-condition types and
              values at the western, eastern, southern, and northern
              boundaries.
    rhs     : Right-hand-side vector.
    K       : Sparse coefficient matrix for the linear system.
    Num     : Structure or tuple containing the global numbering of the
              centroid nodes.

# Keyword Arguments

    Q       : Volumetric heat-production rate (default: `zeros(NC...)`).
    ρ       : Density (default: `3300.0`).
    cp      : Specific heat capacity (default: `1200.0`).
    Qₛ      : Volumetric shear-heating rate (default: `zeros(NC...)`).

# Notes

Backward Euler is first-order accurate in time and unconditionally stable for
the linear diffusion equation. The spatial discretization is second-order
accurate.

The diffusion operator is represented by a five-point finite-difference
stencil. Boundary conditions modify both the matrix coefficients and the
right-hand-side vector for the equations adjacent to the domain boundaries.

The coefficient matrix and right-hand-side vector are assembled inside the
routine. The temperature at the new time level is obtained from

    K * Tⁿ⁺¹ = rhs

and written to `D.T`. The interior of `D.T_ex` is then updated with the new
temperature field.

"""
function BackwardEuler2Dc!(D, κ, Δx, Δy, Δt, NC, BC, rhs, K, Num; 
                        Q = zeros(NC...), ρ = 3300.0, cp = 1200.0, Qₛ = zeros(NC...) )
# dT/dt = kappa*d^2T_ij/dx_i^2 + Q_ij/ρ/cp + Qₛ[i,j]
# ----------------------------------------------------------------------- #
# Define coefficients ---
a   =   κ / Δx^2
b   =   κ / Δy^2
c   =   1 / Δt

rhs  .= reshape(D.T,NC.x*NC.y).*c .+ 
            reshape(Q,NC.x*NC.y)./ρ./cp + reshape(Qₛ,NC.x*NC.y)./ρ./cp

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
    CNA2Dc!(
        D, κ, Δx, Δy, Δt, NC, BC, rhs, K1, K2, Num;
        Q=zeros(NC...), ρ=3300.0, cp=1200.0, Qₛ=zeros(NC...)
    )

Solves the two-dimensional transient heat equation using the Crank–Nicolson
finite-difference scheme with constant thermal diffusivity.

Temperature is defined at the cell centroids, and the diffusion operator is
discretized using second-order central finite differences. Dirichlet and
Neumann boundary conditions are incorporated directly into the coefficient
matrices and right-hand-side vector at the western, eastern, southern, and
northern boundaries.

Optional volumetric heat-production and shear-heating terms may be included.
The resulting linear system is solved using a direct left-matrix division.

# Arguments

    D       : Structure or tuple containing:
              `T`, the temperature field defined at the cell centroids, and
              `T_ex`, the extended temperature field including ghost nodes.
              Both arrays are updated in place.
    κ       : Thermal diffusivity.
    Δx      : Horizontal grid spacing.
    Δy      : Vertical grid spacing.
    Δt      : Time step.
    NC      : Structure or tuple containing the number of centroid nodes
              in the x- and y-directions.
    BC      : Structure or tuple defining the boundary-condition types and
              values at the western, eastern, southern, and northern
              boundaries.
    rhs     : Right-hand-side vector.
    K1      : Sparse coefficient matrix associated with the unknown
              temperature field at the new time level.
    K2      : Sparse coefficient matrix associated with the known
              temperature field at the previous time level.
    Num     : Structure or tuple containing the global numbering of the
              centroid nodes.

# Keyword Arguments

    Q       : Volumetric heat-production rate (default: `zeros(NC...)`).
    ρ       : Density (default: `3300.0`).
    cp      : Specific heat capacity (default: `1200.0`).
    Qₛ      : Volumetric shear-heating rate (default: `zeros(NC...)`).

# Notes

The Crank–Nicolson scheme is second-order accurate in both space and time.
For linear diffusion problems it is unconditionally stable, although large
time steps may introduce weak temporal oscillations.

The diffusion operator is represented by a five-point finite-difference
stencil. Boundary conditions modify both the coefficient matrices and the
right-hand-side vector for the equations adjacent to the domain boundaries.

The matrices `K1` and `K2` represent the implicit and explicit contributions
of the generalized time discretization, respectively. The right-hand-side
vector is assembled from the previous temperature field and the volumetric
heat-source terms before solving

    K1 * Tⁿ⁺¹ = rhs

The computed temperature field is written to `D.T`, and the interior of
`D.T_ex` is updated accordingly.
"""
function CNA2Dc!(D, κ, Δx, Δy, Δt, NC, BC, rhs, K1, K2, Num; 
                Q = zeros(NC...), ρ = 3300.0, cp = 1200.0, Qₛ = zeros(NC...) )
# dT/dt = kappa*d^2T_ij/dx_i^2 + Q_ij/ρ/cp + Qₛ[i,j]/ρ/cp
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
                reshape(Q,NC.x*NC.y)./ρ./cp + reshape(Qₛ,NC.x*NC.y)./ρ./cp
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
                        2*b*BC.val.S[i]*Δy * NeuS + 
                        2*b*BC.val.N[i]*Δy * NeuN
    end
end
# ------------------------------------------------------------------- #    
# Compute new temperature ------------------------------------------- #
D.T[:]      .=    K1 \ rhs[:]
D.T_ex[2:end-1,2:end-1]     .=    D.T
# ------------------------------------------------------------------- #
end

"""
    ADI2Dc!(
        T, κ, Δx, Δy, Δt, NC, BC;
        Q=zeros(NC...), ρ=3300.0, cp=1200.0
    )

Solves the two-dimensional transient heat equation using a second-order
Alternating Direction Implicit (ADI) finite-difference scheme with constant
thermal diffusivity.

Temperature is defined at the cell centroids, and the diffusion operator is
discretized using second-order central finite differences. Dirichlet and
Neumann boundary conditions are incorporated directly into the linear systems
at the western, eastern, southern, and northern boundaries.

Optional volumetric heat production may be included.

# Arguments

    T       : Structure or tuple containing:
              `T`, the temperature field defined at the cell centroids, and
              `T_ex`, the extended temperature field including ghost nodes.
              Both arrays are updated in place.
    κ       : Thermal diffusivity.
    Δx      : Horizontal grid spacing.
    Δy      : Vertical grid spacing.
    Δt      : Time step.
    NC      : Structure or tuple containing the number of centroid nodes
              in the x- and y-directions.
    BC      : Structure or tuple defining the boundary-condition types and
              values at the western, eastern, southern, and northern
              boundaries.

# Keyword Arguments

    Q       : Volumetric heat-production rate (default: `zeros(NC...)`).
    ρ       : Density (default: `3300.0`).
    cp      : Specific heat capacity (default: `1200.0`).

# Notes

The ADI method advances the solution over one time step by performing two
successive implicit half-steps:

1. An implicit solve in the vertical direction while treating the horizontal
   diffusion term explicitly.
2. An implicit solve in the horizontal direction while treating the vertical
   diffusion term explicitly.

The symmetric combination of the two directional half-steps yields
second-order accuracy in time. The spatial diffusion operators are
discretized using second-order central finite differences, making the method
second-order accurate in both space and time.

For the linear two-dimensional diffusion equation, the ADI scheme is
unconditionally stable. Each half-step requires the solution of a system that
is implicit in only one spatial direction, rather than a fully coupled
two-dimensional system.

The temperature field is updated in place, and the interior of `T.T_ex` is
synchronized with the new temperature field after the second half-step.
"""
function ADI2Dc!(T, κ, Δx, Δy, Δt, NC, BC; 
                Q = zeros(NC...), ρ = 3300.0, cp = 1200.0 )
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
                    reshape(Q',(NC.y*NC.x,1))./ρ./cp

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
                            a*BC.val.W[j] * Δx * NeuW + 
                            a*BC.val.E[j] * Δx * NeuE - 
                            b*BC.val.S[i] * Δy * NeuS + 
                            b*BC.val.N[i] * Δy * NeuN    
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
                    reshape(Q,(NC.y*NC.x,1))./ρ./cp
    
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
                            a*BC.val.W[j] * Δx * NeuW + 
                            a*BC.val.E[j] * Δx * NeuE - 
                            b*BC.val.S[i] * Δy * NeuS + 
                            b*BC.val.N[i] * Δy * NeuN    
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
    ComputeResiduals2D!(
        R, T, T_ex, T0, T_ex0, Q, ∂T, q, ρ, Cp, k, BC, Δ, Δt;
        C=0.0, Qₛ=0.0
    )

Computes the residual of the conservative two-dimensional transient heat
equation with variable thermal conductivity and optional volumetric heat
sources.

Temperature is defined at the cell centroids, while heat fluxes are evaluated
at the cell faces using Fourier's law. The divergence of the heat flux is then
computed using second-order finite differences, providing a conservative
discretization of the diffusion operator.

Ghost nodes are used to impose Dirichlet and Neumann boundary conditions at
the western, eastern, southern, and northern boundaries.

The temporal discretization is controlled by `C` and can represent Backward
Euler, Crank–Nicolson, or Forward Euler time integration. The residual is
intended for defect-correction iterations and is used together with the
coefficient matrix assembled by `AssembleMatrix2D`.

# Arguments

    R       : Residual field defined at the cell centroids.
    T       : Temperature field at the current iteration or new time level.
    T_ex    : Extended current-temperature field including ghost nodes.
    T0      : Temperature field at the previous time level.
    T_ex0   : Extended previous-temperature field including ghost nodes.
    Q       : Volumetric heat-production rate.
    ∂T      : Structure or tuple containing the temperature gradients
              `∂x` and `∂y`.
    q       : Structure or tuple containing the heat fluxes
              (`x`, `y`, `x0`, and `y0`).
    ρ       : Density.
    Cp      : Specific heat capacity.
    k       : Thermal conductivity defined on the cell faces.
    BC      : Structure or tuple defining the boundary-condition types and
              values at the western, eastern, southern, and northern
              boundaries.
    Δ       : Structure or tuple containing the horizontal and vertical grid
              spacings `x` and `y`.
    Δt      : Time step.

# Keyword Arguments

    C       : Temporal weighting parameter (default: `0.0`):
              `0.0` for Backward Euler,
              `0.5` for Crank–Nicolson,
              `1.0` for Forward Euler.
    Qₛ      : Volumetric shear-heating rate (default: `0.0`).

# Notes

The residual is evaluated from

- the transient heat-storage term,
- the divergence of the conductive heat flux,
- the volumetric heat-production term, and
- the volumetric shear-heating term.

The heat flux is computed from Fourier's law,

    q = -k ∇T,

using thermal conductivities defined at the cell faces. The divergence of the
heat flux is then evaluated at the cell centroids, resulting in a conservative
finite-difference discretization that naturally accommodates spatially varying
thermal conductivity.

Backward Euler evaluates the conductive fluxes entirely at the new time level,
Crank–Nicolson averages the previous and current time levels, and Forward
Euler evaluates the conductive fluxes at the previous time level.
    
"""
function ComputeResiduals2D!(R, T, T_ex, T0, T_ex0, Q, ∂T, q, ρ, cp, k, BC, Δ, Δt;C=0,Qₛ=0.0)
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
            @. ∂T.∂x0 = (T_ex0[2:end,2:end-1] - T_ex0[1:end-1,2:end-1])/Δ.x
            @. ∂T.∂y0 = (T_ex0[2:end-1,2:end] - T_ex0[2:end-1,1:end-1])/Δ.y
            @. q.x0  = -k.x * ∂T.∂x0
            @. q.y0  = -k.y * ∂T.∂y0
        end
    else
        @. T_ex0[2:end-1,2:end-1] = T0 
        @. T_ex0[  1,2:end-1] = (BC.type.W==:Dirichlet) * (2*BC.val.W - T_ex0[    2,2:end-1]) + (BC.type.W==:Neumann) * (T_ex0[    2,2:end-1] - Δ.x/k.x[  1,:]*BC.val.W)
        @. T_ex0[end,2:end-1] = (BC.type.E==:Dirichlet) * (2*BC.val.E - T_ex0[end-1,2:end-1]) + (BC.type.E==:Neumann) * (T_ex0[end-1,2:end-1] + Δ.x/k.x[end,:]*BC.val.E)
        @. T_ex0[2:end-1,  1] = (BC.type.S==:Dirichlet) * (2*BC.val.S - T_ex0[2:end-1,    2]) + (BC.type.S==:Neumann) * (T_ex0[2:end-1,    2] - Δ.y/k.y[:,  1]*BC.val.S)
        @. T_ex0[2:end-1,end] = (BC.type.N==:Dirichlet) * (2*BC.val.N - T_ex0[2:end-1,end-1]) + (BC.type.N==:Neumann) * (T_ex0[2:end-1,end-1] + Δ.y/k.y[:,end]*BC.val.N)
        @. ∂T.∂x0 = (T_ex0[2:end,2:end-1] - T_ex0[1:end-1,2:end-1])/Δ.x
        @. ∂T.∂y0 = (T_ex0[2:end-1,2:end] - T_ex0[2:end-1,1:end-1])/Δ.y
        @. q.x0  = -k.x * ∂T.∂x0
        @. q.y0  = -k.y * ∂T.∂y0
        @. q.x   = q.x0
        @. q.y   = q.y0
    end
    # @. R     = (T - T0)/Δt - ((1-C)*(∂2T.∂x2 + ∂2T.∂y2) + C*(∂2T.∂x20 + ∂2T.∂y20)) - Q
    @. R     = ρ*cp*(T - T0)/Δt + 
                    (1-C)*((q.x[2:end,:] - q.x[1:end-1,:])/Δ.x + (q.y[:,2:end] - q.y[:,1:end-1])/Δ.y) + 
                    C*((q.x0[2:end,:] - q.x0[1:end-1,:])/Δ.x + (q.y0[:,2:end] - q.y0[:,1:end-1])/Δ.y) - 
                    Q - Qₛ
end

"""
    AssembleMatrix2D(ρ, cp, k, BC, Num, nc, Δ, Δt; C=0.0)

Assembles the coefficient matrix for the conservative two-dimensional
transient heat equation with spatially variable thermal properties.

Temperature is defined at the cell centroids, while thermal conductivity is
defined at the cell faces. The conductive term is discretized in flux form
using second-order finite differences, resulting in a conservative
discretization of the variable-conductivity diffusion operator.

Dirichlet and Neumann boundary conditions are incorporated directly into the
matrix coefficients at the western, eastern, southern, and northern
boundaries.

# Arguments

    ρ       : Density defined at the cell centroids.
    cp      : Specific heat capacity defined at the cell centroids.
    k       : Structure or tuple containing the face-centered thermal
              conductivities `x` and `y`.
    BC      : Structure or tuple defining the boundary-condition types and
              values at the western, eastern, southern, and northern
              boundaries.
    Num     : Structure or tuple containing the global numbering of the
              centroid nodes.
    nc      : Structure or tuple containing the number of centroid nodes
              in the x- and y-directions.
    Δ       : Structure or tuple containing the horizontal and vertical grid
              spacings `x` and `y`.
    Δt      : Time step.

# Keyword Arguments

    C       : Temporal weighting parameter (default: `0.0`):
              `0.0` for Backward Euler,
              `0.5` for Crank–Nicolson,
              `1.0` for Forward Euler.

# Notes

The matrix corresponds to the current-time contribution of the generalized
time-discretization scheme used by `ComputeResiduals2D!`.

The transient storage term is weighted by the local volumetric heat capacity

    ρ * cp,

while the conductive coefficients are determined from the thermal
conductivities at the western, eastern, southern, and northern cell faces.

The spatial operator is represented by a five-point finite-difference
stencil. For an interior cell, the corresponding matrix row couples the
central temperature to its four neighboring centroid temperatures. This
results in a sparse matrix with five non-zero diagonals for the adopted global
numbering.

For `C = 1`, the implicit conductive contribution vanishes, and the matrix
contains only the transient storage term. Volumetric heat sources do not enter
the coefficient matrix and are instead included in the residual or
right-hand-side vector.

The matrix is assembled internally as an `ExtendableSparseMatrix` and
finalized using `flush!(K)` before being returned.
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
    Poisson2Dc!(D, NC, P, BC, Δ, K, rhs, Num)

Solves the two-dimensional steady-state heat equation (Poisson equation)
assuming constant thermal conductivity and volumetric heat production.

Temperature is defined at the cell centroids, and the diffusion operator is
discretized using second-order central finite differences. Dirichlet and
Neumann boundary conditions are incorporated directly into the coefficient
matrix and right-hand-side vector at the western, eastern, southern, and
northern boundaries.

The resulting linear system is solved using a direct left-matrix division.

# Arguments

    D       : Structure or tuple containing the centroid temperature field
              `T` and the volumetric heat-production rate `Q`.
    NC      : Structure or tuple containing the number of centroid nodes
              in the x- and y-directions.
    P       : Structure or tuple containing the thermal conductivity `k`.
    BC      : Structure or tuple defining the boundary-condition types and
              values at the western, eastern, southern, and northern
              boundaries.
    Δ       : Structure or tuple containing the horizontal and vertical grid
              spacings.
    K       : Sparse coefficient matrix.
    rhs     : Right-hand-side vector.
    Num     : Structure or tuple containing the global numbering of the
              centroid nodes.

# Notes

This routine solves the steady-state heat equation where the thermal conductivity 
is assumed to be constant throughout the domain.

The diffusion operator is discretized using a classical five-point finite
difference stencil, resulting in a sparse matrix with five non-zero diagonals
for the adopted global numbering.

Dirichlet and Neumann boundary conditions are incorporated directly into the
coefficient matrix and right-hand-side vector. After assembly, the steady-state
temperature field is obtained from

KT=rhs,

and written to D.T.
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

    # D.T_ex[2:end-1,2:end-1]     .=    D.T
    
end

# ======================================================================= #
# Steady state solvers, variable thermal parameters ===================== #
# ======================================================================= #
"""
    Poisson2D!(T, Q, kx, ky, Δx, Δy, NC, BC, K, rhs, Num)

Solves the conservative two-dimensional steady-state heat equation with
spatially variable thermal conductivity and optional volumetric heat
production.

Temperature is defined at the cell centroids, while the horizontal and
vertical thermal conductivities are defined at the corresponding cell faces.
The conductive term is discretized in flux form using second-order finite
differences.

Dirichlet and Neumann boundary conditions are incorporated directly into the
coefficient matrix and right-hand-side vector at the western, eastern,
southern, and northern boundaries. The resulting linear system is solved using
a direct left-matrix division.

# Arguments

    T       : Temperature field defined at the cell centroids. The array is
              updated in place with the steady-state solution.
    Q       : Volumetric heat-production rate defined at the cell centroids.
    kx      : Thermal conductivity defined at the vertical cell faces.
    ky      : Thermal conductivity defined at the horizontal cell faces.
    Δx      : Horizontal grid spacing.
    Δy      : Vertical grid spacing.
    NC      : Structure or tuple containing the number of centroid nodes
              in the x- and y-directions.
    BC      : Structure or tuple defining the boundary-condition types and
              values at the western, eastern, southern, and northern
              boundaries.
    K       : Sparse coefficient matrix.
    rhs     : Right-hand-side vector.
    Num     : Structure or tuple containing the global numbering of the
              centroid nodes.

# Notes

This routine solves the steady-state heat equation in conservative form,
where thermal conductivity may vary spatially.

The horizontal and vertical conductive fluxes are evaluated using the
face-centered conductivity fields kx and ky. Their divergence is then
discretized at the cell centroids using a five-point finite-difference stencil.
This formulation conserves conductive heat flux and is suitable for material
interfaces with discontinuous thermal conductivity.

Dirichlet and Neumann boundary conditions modify both the matrix coefficients
and the right-hand-side vector for equations adjacent to the domain
boundaries.

After assembly, the steady-state temperature field is obtained from

KT=rhs,

and written to T.
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
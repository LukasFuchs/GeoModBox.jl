function ForwardEuler_const!(D,NC,BC,P,Δ,T)
    # Function to solve 2D heat diffusion equation using the explicit finite
    # difference scheme
    # Q - Waermeproduktionsrate pro Volumen [W/m^3]
    # ------------------------------------------------------------------- #
    
    sx      = P.κ[1] * T.Δ[1] / Δ.x[1]^2
    sz      = P.κ[1] * T.Δ[1] / Δ.y[1]^2

    D.T_ex[2:end-1,2:end-1]     .=  D.T

    # Temperature at the ghost nodes ------------------------------------ #
    # West boundary ---
    D.T_ex[1,2:end-1]   .= (BC.type.W==:Dirichlet) .* (2 .* BC.val.W .- D.T_ex[2,2:end-1])
    # East boundary ---
    D.T_ex[end,2:end-1] .= (BC.type.E==:Dirichlet) .* (2 .* BC.val.E .- D.T_ex[end-1,2:end-1])
    # South boundary --- 
    D.T_ex[2:end-1,1]   .= (BC.type.S==:Dirichlet) .* (2 .* BC.val.S .- D.T_ex[2:end-1,2])
    # Northern boundary ---
    D.T_ex[2:end-1,end] .= (BC.type.N==:Dirichlet) .* (2 .* BC.val.N .- D.T_ex[2:end-1,end-1])
    # ------------------------------------------------------------------- #
    
    # Loop over internal nodes ------------------------------------------ #
    for i = 1:NC.x, j = 1:NC.y
        i1 = i+1
        j1 = j+1
        D.T[i,j] = D.T_ex[i1,j1] + 
            sx * (D.T_ex[i1-1,j1] - 2 * D.T_ex[i1,j1] + D.T_ex[i1+1,j1]) + 
            sz * (D.T_ex[i1,j1-1] - 2 * D.T_ex[i1,j1] + D.T_ex[i1,j1+1]) + 
            D.Q[i,j] * T.Δ[1] / P.ρ / P.cp
    end
    # ------------------------------------------------------------------- #
    
    end
using ExtendableSparse

function Assemblyc(NC, NV, Δ, η, BC, Num)

    # Linear system of equation ---
    ndof    =   maximum(Num.Pt)
    K       =   ExtendableSparseMatrix(ndof,ndof)
    dx,dy   =   Δ.x, Δ.y

    # x momentum equation ----------------------------------------------- #
    for i = 1:NV.x, j = 1:NC.y
        # Equation Number ---
        ii  =   Num.Vx[i,j] 
        if i == 1 || i == NV.x
            # East and West boundary ---
            # Free Slip && No Slip: vₓ = 0 
            K[ii,ii]    =   1.0
        else
            # Stencil, vₓ ---
            iS  =   ii - NV.x
            iW  =   ii - 1 
            iC  =   ii 
            iE  =   ii + 1
            iN  =   ii + NV.x
            # Pressure ---
            iPE =   Num.Pt[i,j]
            iPW =   Num.Pt[i-1,j]
            # ---
            inS     =   j==1    ? false  : true  
            FSS     =   (j==1    && BC.type.S==:freeslip) ? 1. : 0.
            NSS     =   (j==1    && (BC.type.S==:noslip||BC.type.S==:const)) ? 1. : 0.
            inN     =   j==NC.y ? false  : true   
            FSN     =   (j==NC.y && BC.type.N==:freeslip) ? 1. : 0.
            NSN     =   (j==NC.y && (BC.type.N==:noslip||BC.type.N==:const)) ? 1. : 0.
            if inS K[ii,iS] =   η / dy^2 end
            K[ii,iW]    =   η / dx^2
            K[ii,iC]    =   - 2.0 * η / dx^2 - (2.0 - FSS - FSN + NSS + NSN) * η / dy^2
            K[ii,iE]    =   η / dx^2
            if inN K[ii,iN] =   η / dy^2 end
            K[ii,iPW]   = 1 / dx
            K[ii,iPE]   = -1 / dx
        end
    end
    # ------------------------------------------------------------------- #
    # y momentum equation ----------------------------------------------- #
    for i = 1:NC.x, j = 1:NV.y
        # Equation Number ---
        ii  =   Num.Vy[i,j] 
        if j == 1 || j == NV.y
            # East and West boundary ---
            # Free Slip && No Slip: vy = 0 
            K[ii,ii]    =   1.0
        else
            # Stencil, vy ---
            iS  =   ii - NC.x
            iW  =   ii - 1 
            iC  =   ii 
            iE  =   ii + 1
            iN  =   ii + NC.x
            # Pressure ---
            iPS =   Num.Pt[i,j-1]   
            iPN =   Num.Pt[i,j]
            # ---
            inW     =   i==1    ? false  : true  
            FSW     =   (i==1    && BC.type.W==:freeslip) ? 1. : 0.
            NSW     =   (i==1    && (BC.type.W==:noslip||BC.type.W==:const)) ? 1. : 0.
            inE     =   i==NC.x ? false  : true   
            FSE     =   (i==NC.x && BC.type.E==:freeslip) ? 1. : 0.
            NSE     =   (i==NC.x && (BC.type.E==:noslip||BC.type.E==:const)) ? 1. : 0.            
            K[ii,iS]    =   η / dy^2
            if inW K[ii,iW] =  η / dx^2 end
            K[ii,iC]    =   - (2.0 - FSW - FSE + NSW + NSE) * η / dx^2 - 2.0 * η / dy^2
            if inE K[ii,iE] = η / dx^2 end
            K[ii,iN]    =   η / dy^2
            K[ii,iPS]   =   1 / dy
            K[ii,iPN]   =   - 1 / dy
        end
    end
    # ------------------------------------------------------------------- #
    # Continuity equation ----------------------------------------------- #
    for i = 1:NC.x, j = 1:NC.y
        # Equation number ---
        ii  =   Num.Pt[i,j]
        # Stencil ---
        iW  =   Num.Vx[i,j]
        iE  =   Num.Vx[i+1,j]
        iS  =   Num.Vy[i,j]
        iN  =   Num.Vy[i,j+1]
        # Linear system coefficients
        K[ii,iW]    =   -1 / dx
        K[ii,iE]    =   1 / dx
        K[ii,iS]    =   -1 / dy 
        K[ii,iN]    =   1 / dy
    end
    return flush!(K)
end

function updaterhsc(NC, NV, Δ, η, ρ, g, BC, Num)

    rhs     =   zeros(maximum(Num.Pt))  #   Right-hand Side

    # x momentum equation ----------------------------------------------- #
    for i = 1:NV.x, j = 1:NC.y
        # Equation Number ---
        ii  =   Num.Vx[i,j] 
        rhs[ii]     =   0.0
        if i == 1 || i == NV.x
            # East and West boundary ---
            # Free Slip && No Slip: vₓ = 0 
            # rhs[ii]     =   0.0
            CW  =   (i==1 && BC.type.W==:const) ? 1. : 0.
            CE  =   (i==NV.x && BC.type.E==:const) ? 1. : 0. 
            rhs[ii]  += CW * BC.val.vxW[j] + CE * BC.val.vxE[j]
        else            
            # ---
            NSS     =   (j==1    && (BC.type.S==:noslip||BC.type.S==:const)) ? 1. : 0.
            NSN     =   (j==NC.y && (BC.type.N==:noslip||BC.type.N==:const)) ? 1. : 0.            
            # ---
            rhs[ii] +=  -2.0 * η * BC.val.S[i] / Δ.y^2 * NSS -
                        2.0 * η * BC.val.N[i] / Δ.y^2 * NSN
        end
    end
    # y momentum equation ----------------------------------------------- #
    for i = 1:NC.x, j = 1:NV.y
        # Equation Number ---
        ii  =   Num.Vy[i,j] 
        rhs[ii]     =   0.0
        if j == 1 || j == NV.y
            # North and South boundary ---
            # Free Slip && No Slip: vy = 0 
            # rhs[ii]     =   0.0
            CS  =   (j==1 && BC.type.S==:const) ? 1. : 0.
            CN  =   (j==NV.y && BC.type.N==:const) ? 1. : 0.
            rhs[ii]     += CN * BC.val.vyN[i] + CS * BC.val.vyS[i]
        else            
            # ---
            NSW     =   (i==1    && (BC.type.W==:noslip||BC.type.W==:const)) ? 1. : 0.
            NSE     =   (i==NC.x && (BC.type.E==:noslip||BC.type.E==:const)) ? 1. : 0.            
            # ---
            rhs[ii] +=  - g * ((ρ[i,j] + ρ[i,j-1]) / 2.0) - 
                        2.0 * η * BC.val.W[j] / Δ.x^2 * NSW -
                        2.0 * η * BC.val.E[j] / Δ.x^2 * NSE
        end
    end
    return rhs
end

function Residuals2Dc!(D,BC,ε,τ,divV,Δ,η,g,Fm,FPt)
    @. D.vx[:,1]    = (BC.type.S==:freeslip)*D.vx[:,2]     + (BC.type.S==:noslip||BC.type.S==:const)*(2*BC.val.S - D.vx[:,2])
    @. D.vx[:,end]  = (BC.type.N==:freeslip)*D.vx[:,end-1] + (BC.type.N==:noslip||BC.type.N==:const)*(2*BC.val.N - D.vx[:,end-1])
    @. D.vy[1,:]    = (BC.type.W==:freeslip)*D.vy[2,:]     + (BC.type.W==:noslip||BC.type.W==:const)*(2*BC.val.W - D.vy[2,:])
    @. D.vy[end,:]  = (BC.type.E==:freeslip)*D.vy[end-1,:] + (BC.type.E==:noslip||BC.type.E==:const)*(2*BC.val.E - D.vy[end-1,:])
    # @. D.vx[:,1]    = (BC.type.S==:freeslip)*D.vx[:,2]     + (BC.type.S==:noslip)*(2*BC.val.S - D.vx[:,2])
    # @. D.vx[:,end]  = (BC.type.N==:freeslip)*D.vx[:,end-1] + (BC.type.N==:noslip)*(2*BC.val.N - D.vx[:,end-1])
    # @. D.vy[1,:]    = (BC.type.W==:freeslip)*D.vy[2,:]     + (BC.type.W==:noslip)*(2*BC.val.W - D.vy[2,:])
    # @. D.vy[end,:]  = (BC.type.E==:freeslip)*D.vy[end-1,:] + (BC.type.E==:noslip)*(2*BC.val.E - D.vy[end-1,:])
    @. divV =   (D.vx[2:end,2:end-1] - D.vx[1:end-1,2:end-1])/Δ.x + (D.vy[2:end-1,2:end] - D.vy[2:end-1,1:end-1])/Δ.y
    @. ε.xx =   (D.vx[2:end,2:end-1] - D.vx[1:end-1,2:end-1])/Δ.x - 1.0/3.0*divV
    @. ε.yy =   (D.vy[2:end-1,2:end] - D.vy[2:end-1,1:end-1])/Δ.y - 1.0/3.0*divV
    @. ε.xy =   0.5*( (D.vx[:,2:end] - D.vx[:,1:end-1])/Δ.y + (D.vy[2:end,:] - D.vy[1:end-1,:])/Δ.x ) 
    @. τ.xx =   2.0 * η * ε.xx
    @. τ.yy =   2.0 * η * ε.yy
    @. τ.xy =   2.0 * η * ε.xy
    @. Fm.x[2:end-1,:] = ((τ.xx[2:end,:]-τ.xx[1:end-1,:])/Δ.x + (τ.xy[2:end-1,2:end]-τ.xy[2:end-1,1:end-1])/Δ.y - (D.Pt[2:end,:]-D.Pt[1:end-1,:])/Δ.x)
    @. Fm.y[:,2:end-1] = ((τ.yy[:,2:end]-τ.yy[:,1:end-1])/Δ.y + (τ.xy[2:end,2:end-1]-τ.xy[1:end-1,2:end-1])/Δ.x - (D.Pt[:,2:end]-D.Pt[:,1:end-1])/Δ.y - g * ((D.ρ[:,2:end] + D.ρ[:,1:end-1]) / 2.0))
    @. FPt             = divV
end

function Assembly(NC, NV, Δ, ηc, ηv, BC, Num)

    # Linear system of equation ---
    ndof    =   maximum(Num.Pt)
    K       =   ExtendableSparseMatrix(ndof,ndof)
    dx,dy   =   Δ.x, Δ.y

    # x momentum equation ----------------------------------------------- #
    for i = 1:NV.x, j = 1:NC.y
        # Equation Number ---
        ii  =   Num.Vx[i,j] 
        if i == 1 || i == NV.x
            # East and West boundary ---
            # Free Slip && No Slip: vₓ = 0 
            K[ii,ii]    =   1.0
        else
            # Stencil, vₓ ---
            iS  =   ii - NV.x            
            iW  =   ii - 1 
            iC  =   ii 
            iE  =   ii + 1            
            iN  =   ii + NV.x
            # --- 
            iSW =   Num.Vy[i-1,j]
            iSE =   iSW + 1
            iNW =   iSW + NC.x
            iNE =   iSW + NC.x + 1
            # iNW =   iSW + NC.y
            # iNE =   iSW + NC.y + 1
            # Pressure ---
            iPC =   Num.Pt[i,j]
            iPW =   Num.Pt[i-1,j]
            # ---
            inS     =   j==1    ? false  : true  
            FSS     =   (j==1    && BC.type.S==:freeslip) ? 1. : 0.
            NSS     =   (j==1    && (BC.type.S==:noslip||BC.type.S==:const)) ? 1. : 0.
            inN     =   j==NC.y ? false  : true   
            FSN     =   (j==NC.y && BC.type.N==:freeslip) ? 1. : 0.
            NSN     =   (j==NC.y && (BC.type.N==:noslip||BC.type.N==:const)) ? 1. : 0.
            if inS K[ii,iS] =   ηv[i,j] / dy^2   end
            K[ii,iSW]   =   ηv[i,j] / dx / dy 
            K[ii,iSE]   =   - ηv[i,j] / dx / dy
            K[ii,iW]    =   2* ηc[i-1,j] / dx^2
            K[ii,iC]    =   -2*ηc[i,j] / dx^2 - 2*ηc[i-1,j] /dx^2 - (1.0-FSN+NSN)*ηv[i,j+1]/dy^2 - (1.0-FSS+NSS)*ηv[i,j]/dy^2
            K[ii,iE]    =   2 * ηc[i,j] / dx^2
            K[ii,iNW]   =   - ηv[i,j+1] / dx / dy
            K[ii,iNE]   =   ηv[i,j+1] / dx / dy
            if inN K[ii,iN] =   ηv[i,j+1] / dy^2  end
            K[ii,iPW]   =   1 / dx
            K[ii,iPC]   =   -1 / dx
        end
    end
    # ------------------------------------------------------------------- #
    # y momentum equation ----------------------------------------------- #
    for i = 1:NC.x, j = 1:NV.y
        # Equation Number ---
        ii  =   Num.Vy[i,j] 
        if j == 1 || j == NV.y
            # North and South boundary ---
            # Free Slip && No Slip: vₓ = 0 
            K[ii,ii]    =   1.0
        else
            # Stencil, vₓ ---
            iS  =   ii - NC.x
            iW  =   ii - 1 
            iC  =   ii 
            iE  =   ii + 1
            iN  =   ii + NC.x
            # --- 
            iSW =   Num.Vx[i,j-1]
            iSE =   iSW + 1
            iNW =   iSW + NV.x
            iNE =   iSW + NV.x + 1
            # Pressure ---
            iPS =   Num.Pt[i,j-1]   
            iPC =   Num.Pt[i,j]
            # ---
            inW     =   i==1    ? false  : true  
            FSW     =   (i==1    && BC.type.W==:freeslip) ? 1. : 0.
            NSW     =   (i==1    && (BC.type.W==:noslip||BC.type.W==:const)) ? 1. : 0.
            inE     =   i==NC.x ? false  : true   
            FSE     =   (i==NC.x && BC.type.E==:freeslip) ? 1. : 0.
            NSE     =   (i==NC.x && (BC.type.E==:noslip||BC.type.E==:const)) ? 1. : 0.            
            K[ii,iS]    =   2 * ηc[i,j-1] / dy^2
            K[ii,iSW]   =   ηv[i,j] / dx / dy
            K[ii,iSE]   =   - ηv[i+1,j] / dx / dy
            if inW K[ii,iW] =  ηv[i,j] / dx^2 end
            K[ii,iC]    =   -2.0*ηc[i,j]/dy^2 - 2.0*ηc[i,j-1]/dy^2 - (1.0-FSE+NSE) * ηv[i+1,j] / dx^2 - (1.0-FSW+NSW) * ηv[i,j] / dx^2
            if inE K[ii,iE] = ηv[i+1,j] / dx^2 end
            K[ii,iNW]   =   - ηv[i,j] / dx / dy
            K[ii,iNE]   =   ηv[i+1,j] / dx / dy
            K[ii,iN]    =   2 * ηc[i,j] / dy^2
            K[ii,iPS]   =   1 / dy
            K[ii,iPC]   =   - 1 / dy
        end
    end
    # ------------------------------------------------------------------- #
    # Continuity equation ----------------------------------------------- #
    for i = 1:NC.x, j = 1:NC.y
        # Equation number ---
        ii  =   Num.Pt[i,j]
        # Stencil ---
        iW  =   Num.Vx[i,j]
        iE  =   Num.Vx[i+1,j]
        iS  =   Num.Vy[i,j]
        iN  =   Num.Vy[i,j+1]
        # Linear system coefficients
        K[ii,iW]    =   -1 / dx
        K[ii,iE]    =   1 / dx
        K[ii,iS]    =   -1 / dy 
        K[ii,iN]    =   1 / dy
    end
    return flush!(K)
end

function updaterhs(NC, NV, Δ, ηc, ηv, ρ, g, BC, Num)

    rhs     =   zeros(maximum(Num.Pt))  #   Right-hand Side

    # x momentum equation ----------------------------------------------- #
    for i = 1:NV.x, j = 1:NC.y
        # Equation Number ---
        ii  =   Num.Vx[i,j] 
        rhs[ii]     =   0.0
        if i == 1 || i == NV.x
            # East and West boundary ---
            # Free Slip && No Slip: vₓ = 0 
            # rhs[ii]     =   0.0
            CW  =   (i==1 && BC.type.W==:const) ? 1. : 0.
            CE  =   (i==NV.x && BC.type.E==:const) ? 1. : 0. 
            rhs[ii]  += CW * BC.val.vxW[j] + CE * BC.val.vxE[j]
        else            
            # ---
            NSS     =   (j==1    && (BC.type.S==:noslip||BC.type.S==:const)) ? 1. : 0.
            NSN     =   (j==NC.y && (BC.type.N==:noslip||BC.type.N==:const)) ? 1. : 0.            
            # ---
            rhs[ii] +=  -2.0 * ηv[i,j] * BC.val.S[i] / Δ.y^2 * NSS -
                        2.0 * ηv[i,j+1] * BC.val.N[i] / Δ.y^2 * NSN
        end
    end
    # y momentum equation ----------------------------------------------- #
    for i = 1:NC.x, j = 1:NV.y
        # Equation Number ---
        ii  =   Num.Vy[i,j] 
        rhs[ii]     =   0.0
        if j == 1 || j == NV.y
            # North and South boundary ---
            # Free Slip && No Slip: vy = 0 
            # rhs[ii]     =   0.0
            CS  =   (j==1 && BC.type.S==:const) ? 1. : 0.
            CN  =   (j==NV.y && BC.type.N==:const) ? 1. : 0.
            rhs[ii]     += CN * BC.val.vyN[i] + CS * BC.val.vyS[i]
        else            
            # ---
            NSW     =   (i==1    && (BC.type.W==:noslip||BC.type.W==:const)) ? 1. : 0.
            NSE     =   (i==NC.x && (BC.type.E==:noslip||BC.type.E==:const)) ? 1. : 0.            
            # ---
            rhs[ii] +=  g * ((ρ[i,j] + ρ[i,j-1]) / 2.0) - 
                        2.0 * ηv[i,j] * BC.val.W[j] / Δ.x^2 * NSW -
                        2.0 * ηv[i+1,j] * BC.val.E[j] / Δ.x^2 * NSE
        end
    end
    return rhs
end

function Residuals2D!(D,BC,ε,τ,divV,Δ,ηc,ηv,g,Fm,FPt)
    @. D.vx[:,1]    = (BC.type.S==:freeslip)*D.vx[:,2]     + (BC.type.S==:noslip||BC.type.S==:const)*(2*BC.val.S - D.vx[:,2])
    @. D.vx[:,end]  = (BC.type.N==:freeslip)*D.vx[:,end-1] + (BC.type.N==:noslip||BC.type.N==:const)*(2*BC.val.N - D.vx[:,end-1])
    @. D.vy[1,:]    = (BC.type.W==:freeslip)*D.vy[2,:]     + (BC.type.W==:noslip||BC.type.W==:const)*(2*BC.val.W - D.vy[2,:])
    @. D.vy[end,:]  = (BC.type.E==:freeslip)*D.vy[end-1,:] + (BC.type.E==:noslip||BC.type.E==:const)*(2*BC.val.E - D.vy[end-1,:])
    @. divV =   (D.vx[2:end,2:end-1] - D.vx[1:end-1,2:end-1])/Δ.x + (D.vy[2:end-1,2:end] - D.vy[2:end-1,1:end-1])/Δ.y
    @. ε.xx =   (D.vx[2:end,2:end-1] - D.vx[1:end-1,2:end-1])/Δ.x - 1.0/3.0*divV
    @. ε.yy =   (D.vy[2:end-1,2:end] - D.vy[2:end-1,1:end-1])/Δ.y - 1.0/3.0*divV
    @. ε.xy =   0.5*( (D.vx[:,2:end] - D.vx[:,1:end-1])/Δ.y + (D.vy[2:end,:] - D.vy[1:end-1,:])/Δ.x ) 
    @. τ.xx =   2.0 * ηc * ε.xx
    @. τ.yy =   2.0 * ηc * ε.yy
    @. τ.xy =   2.0 * ηv * ε.xy
    @. Fm.x[2:end-1,:] = ((τ.xx[2:end,:]-τ.xx[1:end-1,:])/Δ.x + (τ.xy[2:end-1,2:end]-τ.xy[2:end-1,1:end-1])/Δ.y - (D.Pt[2:end,:]-D.Pt[1:end-1,:])/Δ.x)
    @. Fm.y[:,2:end-1] = ((τ.yy[:,2:end]-τ.yy[:,1:end-1])/Δ.y + (τ.xy[2:end,2:end-1]-τ.xy[1:end-1,2:end-1])/Δ.x - (D.Pt[:,2:end]-D.Pt[:,1:end-1])/Δ.y - g * ((D.ρ[:,2:end] + D.ρ[:,1:end-1]) / 2.0))
    @. FPt             = divV
end
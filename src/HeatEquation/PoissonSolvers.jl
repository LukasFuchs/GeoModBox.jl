using ExtendableSparse

function Poisson!(D,NC,P,BC,Δ,K,rhs,Num)
# Function to solve 2D heat diffusion equation using the explicit finite
# difference scheme
# [Q] = W/m^3
# ----------------------------------------------------------------------- #
    
    # a       =   1.0 / Δ.x[1]^2.0
    # b       =   1.0 / Δ.y[1]^2.0

    #  --------------------------------------------- #
    rhs     .=   - reshape(D.Q, NC.x*NC.y, 1) ./ P.k

    for i=1:NC.x, j=1:NC.y
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

    D.T[:]  .=   K \ rhs[:]
    
end

#function PoissonV!()
#    # Function to solve 2D heat diffusion equation using the explicit finite
#    # difference scheme
#    # [Q] = W/m^3
#    # ----------------------------------------------------------------------- #
#    
#    ## Define Coeffizients for Matrix A and boundary conditions ------------- #
#    
#    
#    # Erstellung der durchlaufenden Indizes ----------------------------- #
#    Number  = zeros(N.nz,N.nx);
#    Number2 = zeros(N.nz,N.nx);
#    num     = 1;
#    for i=1:N.nz
#        for j=1:N.nx
#            l   = j+(i-1)*N.nx;
#            Number(i,j) = num;
#            Number2(i,j) = l;
#            num = num+1;
#        end
#    end
#    
#    # Setup coefficient matrix A ---------------------------------------- #
#    a       = 1/N.dz^2;
#    b       = 1/N.dx^2;
#    # Reshape k into a vector so that the index fits later in the coefficient 
#    # calculations
#    k       = reshape(k',[N.nx*N.nz,1]);
#    
#    # Define diagonals for matrix
#    diag    = zeros(N.nx*N.nz,5);
#    
#    # Inner index ------------------------------------------------------- #
#    iind   = Number(2:(N.nz-1),2:(N.nx-1));
#    
#    diag(iind-N.nx,1)   = a.*(k(iind-N.nx,1)+k(iind,1))./2;
#    diag(iind-1,2)      = b.*(k(iind,1)+k(iind-1,1))./2;
#    diag(iind,3)        = ...
#            -b.*((k(iind+1,1)+k(iind,1))./2+(k(iind,1)+k(iind-1,1))./2) ...
#            - a.*((k(iind+N.nx,1)+k(iind,1))./2+(k(iind-N.nx,1)+k(iind-1,1))./2);
#    diag(iind+1,4)      = b.*(k(iind+1,1)+k(iind,1))./2;
#    diag(iind+N.nx,5)   = a.*(k(iind+N.nx,1)+k(iind,1))./2;
#    
#    # Top and bottom index and boundary conditions ---------------------- #
#    switch B.btbc
#        case {'const','dirichlet'}
#            bind                = Number(1,:);
#            diag(bind,3)        = 1;
#    #     case {'flux','neumann'}
#    #         bind                = Number(1,2:N.nx-1);        
#    #         diag(bind-1,2)      = b;
#    #         diag(bind,3)        = c;
#    #         diag(bind+1,4)      = b;        
#    #         diag(bind+N.nx,5)   = 2*a;   
#    #         
#    #         # Left corner
#    #         blcind              = Number(1,1);        
#    #         diag(blcind,3)      = c;
#    #         diag(blcind+1,4)    = 2*b;        
#    #         diag(blcind+N.nx,5) = 2*a; 
#    #         
#    #         # Right corner
#    #         brcind              = Number(1,N.nx);         
#    #         diag(brcind-1,2)    = 2*b; 
#    #         diag(brcind,3)      = c;   
#    #         diag(blcind+N.nx,5) = 2*a; 
#        otherwise
#            error('Thermal boundary condition not defined!')
#    end
#    switch B.ttbc
#        case {'const','dirichlet'}
#            tind                = Number(N.nz,:);
#            diag(tind,3)        = 1;
#    #     case {'flux','neumann'}
#    #         tind                = Number(N.nz,2:N.nx-1);
#    #         diag(tind-N.nx,1)   = 2*a;
#    #         diag(tind-1,2)      = b;
#    #         diag(tind,3)        = c;
#    #         diag(tind+1,4)      = b;        
#    #         
#    #         # Left corner
#    #         tlcind              = Number(N.nz,1); 
#    #         diag(tlcind-N.nx,1) = 2*a;
#    #         diag(tlcind,3)      = c;
#    #         diag(tlcind+1,4)    = 2*b;
#    #         
#    #         # Right corner
#    #         trcind              = Number(N.nz,N.nx); 
#    #         diag(trcind-N.nx,1) = 2*a;
#    #         diag(trcind-1,2)    = 2*b;
#    #         diag(trcind,3)      = c;        
#        otherwise
#            error('Thermal boundary condition not defined!')
#    end
#    
#    # Left and right index and boundary conditions ---------------------- #
#    lind    = Number(2:(N.nz-1),1);
#    rind    = Number(2:(N.nz-1),N.nx);
#    
#    switch B.ltbc
#        case {'const','dirichlet'}
#            diag(lind,3)    = 1;
#    #     case {'flux','neumann'}
#    #         diag(lind-1,2)      = a;
#    #         diag(lind,3)        = c;
#    #         diag(lind+1,4)      = 2*b;
#    #         diag(lind+N.nx,5)   = a;
#        otherwise
#            error('Thermal boundary condition not defined')
#    end
#    switch B.rtbc
#        case {'const','dirichlet'}
#            diag(rind,3)    = 1;
#    #     case {'flux','neumann'}
#    #         diag(rind-N.nx,1)   = a;
#    #         diag(rind-1,2)      = 2*b;
#    #         diag(rind,3)        = c;
#    #         diag(rind+1,4)      = a;        
#        otherwise
#            error('Thermal boundary condition not defined')
#    end
#    
#    A       = spdiags(diag,[-N.nx,-1,0,1,N.nx],N.nx*N.nz,N.nx*N.nz);
#    
#    ## Solve system of equations -------------------------------------------- #
#    rhs    = -1.*reshape(Q',[N.nx*N.nz,1]);
#    
#    # switch B.btbc
#    #     case {'flux','neumann'}
#    #         #i=1
#    #         rhs(bind)   = rhs(bind) - 2*N.dz*B.bhf*a;        
#    #         rhs(blcind) = rhs(blcind) - 2*N.dx*B.lhf*b - 2*N.dz*B.bhf*a; 
#    #         rhs(brcind) = rhs(brcind) - 2*N.dx*B.rhf*b + 2*N.dz*B.bhf*a; 
#    # end
#    # switch B.ttbc
#    #     case {'flux','neumann'}
#    #         #i=nz
#    #         rhs(tind)   = rhs(tind) - 2*N.dz*B.thf*a;
#    #         rhs(tlcind) = rhs(tlcind) - 2*N.dx*B.lhf*b - 2*N.dz*B.thf*a; 
#    #         rhs(trcind) = rhs(trcind) - 2*N.dx*B.rhf*b - 2*N.dz*B.thf*a; 
#    # end
#    # switch B.ltbc
#    #     case {'flux','neumann'}
#    #         rhs(lind)   = rhs(lind) - 2*N.dx*B.lhf*b;
#    # end
#    # switch B.rtbc
#    #     case {'flux','neumann'}
#    #         rhs(rind)   = rhs(rind) - 2*N.dx*B.rhf*b;
#    # end
#    
#    F      = A\rhs;
#    
#    F      = reshape(F,[N.nx,N.nz])';
#    
#    end
using ExtendableSparse

function ComputeResiduals!(R, T, T_ex, T0, ∂T, q, ρ, Cp, k, BC, Δ, Δt)
    @. T_ex[2:end-1,2:end-1] = T 
    @. T_ex[  1,2:end-1] = (BC.type.W==:Dirichlet) * (2*BC.val.W - T_ex[    2,2:end-1]) + (BC.type.W==:Neumann) * (T_ex[    2,2:end-1] - Δ.x/k.x[  1,:]*BC.val.W)
    @. T_ex[end,2:end-1] = (BC.type.E==:Dirichlet) * (2*BC.val.E - T_ex[end-1,2:end-1]) + (BC.type.E==:Neumann) * (T_ex[end-1,2:end-1] + Δ.x/k.x[end,:]*BC.val.E)
    @. T_ex[2:end-1,  1] = (BC.type.S==:Dirichlet) * (2*BC.val.S - T_ex[2:end-1,    2]) + (BC.type.S==:Neumann) * (T_ex[2:end-1,    2] - Δ.y/k.y[:,  1]*BC.val.S)
    @. T_ex[2:end-1,end] = (BC.type.N==:Dirichlet) * (2*BC.val.N - T_ex[2:end-1,end-1]) + (BC.type.N==:Neumann) * (T_ex[2:end-1,end-1] - Δ.y/k.y[:,end]*BC.val.N)
    @. ∂T.∂x = (T_ex[2:end,2:end-1] - T_ex[1:end-1,2:end-1])/Δ.x
    @. ∂T.∂y = (T_ex[2:end-1,2:end] - T_ex[2:end-1,1:end-1])/Δ.y
    @. q.x   = -k.x * ∂T.∂x
    @. q.y   = -k.y * ∂T.∂y
    @. R     = ρ*Cp*(T - T0)/Δt + (q.x[2:end,:] - q.x[1:end-1,:])/Δ.x + (q.y[:,2:end] - q.y[:,1:end-1])/Δ.y  
end

function AssembleMatrix(rho, cp, k, BC, Num, nc, Δ, Δt)
    # Linear system of equation
    ndof   = maximum(Num.T)
    K      = ExtendableSparseMatrix(ndof, ndof)
    dx, dy = Δ.x, Δ.y
    #############################
    #       Heat equation       #
    #############################
    for i=1:nc.x, j=1:nc.y
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
        kW = k.x[i,j]
        kE = k.x[i+1,j]
        kS = k.y[i,j]
        kN = k.y[i,j+1]
        ρ  = rho[i,j]
        Cp = cp[i,j]
        # Linear system coefficients
        if inS K[ii,iS] = kS .* (DirS + NeuS - 1) ./ dy .^ 2 end
        if inW K[ii,iW] = kW .* (DirW + NeuW - 1) ./ dx .^ 2 end
        K[ii,iC] = Cp .* ρ ./ Δt + (-kN .* (-DirN + NeuN - 1) ./ dy + kS .* (DirS - NeuS + 1) ./ dy) ./ dy + (-kE .* (-DirE + NeuE - 1) ./ dx + kW .* (DirW - NeuW + 1) ./ dx) ./ dx
        if inE K[ii,iE] = -kE .* (-DirE - NeuE + 1) ./ dx .^ 2 end
        if inN K[ii,iN] = -kN .* (-DirN - NeuN + 1) ./ dy .^ 2 end
    end
    return flush!(K)
end

# function Direct!()
# Function to solve 2D heat diffusion equation using the explicit finite
# difference scheme
# assuming constant k, rho, cp
# dT/dt = kappa*d^2T_ij/dx_i^2 + Q_ij/rho/cp
# ----------------------------------------------------------------------- #
#
### Define Coeffizients for Matrix A and boundary conditions ------------- #
## Setup coefficient matrix A -------------------------------------------- #
#sx      =   P.kappa*dt/N.dx^2;
#sz      =   P.kappa*dt/N.dz^2;
#
#if N.beenhere == 0
#    # Erstellung der durchlaufenden Indizes ----------------------------- #
#    Number  =   zeros(N.nz,N.nx);
#    num     =   1;
#    for i=1:N.nz
#        for j=1:N.nx
#            Number(i,j) =   num;
#            num = num+1;
#        end
#    end
#    
#    # Define diagonals for matrix
#    diag    =   zeros(N.nx*N.nz,5);
#    
#    # Inner index ------------------------------------------------------- #
#    iind    =    Number(2:(N.nz-1),2:(N.nx-1));
#    
#    diag(iind-N.nx,1)   =   -sz;
#    diag(iind-1,2)      =   -sx;
#    diag(iind,3)        =   1+2*sx+2*sz;
#    diag(iind+1,4)      =   -sx;
#    diag(iind+N.nx,5)   =   -sz;
#    
#    # Top and bottom index and boundary conditions ---------------------- #
#    switch B.ttbc
#        case {'const','dirichlet'}
#            tind                =   Number(N.nz,:);
#            diag(tind,3)        =   1;
#        case {'flux','neumann'}
#            tind                =   Number(N.nz,2:(N.nx-1));
#            diag(tind-N.nx,1)   =   -2*sz;
#            diag(tind-1,2)      =   -sx;
#            diag(tind,3)        =   1+2*sx+2*sz;
#            diag(tind+1,4)      =   -sx;
#            
#            diag(N.nx*(N.nz-1)+1,3) =   1; # Left Corner grid point
#            diag(N.nx*N.nz,3)       =   1; # Right Corner grid point
#        otherwise
#            error('Thermal boundary condition not defined!')
#    end
#    switch B.btbc
#        case {'const','dirichlet'}
#            bind                =   Number(1,:);
#            diag(bind,3)        =   1;
#        case {'flux','neumann'}
#            bind                =   Number(1,2:(N.nx-1));
#            diag(bind-1,2)      =   -sx;
#            diag(bind,3)        =   1+2*sx+2*sz;
#            diag(bind+1,4)      =   -sx;
#            diag(bind+N.nx,5)   =   -2*sz;
#            
#            diag(1,3)           =   1; # Left Corner grid point
#            diag(N.nx,3)        =   1; # Right Corner grid point
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
#            diag(lind,3)        = 1;
#        case {'flux','neumann'}
#            diag(lind-N.nx,1)   = -sz;
#            diag(lind,3)        = 1+2*sx+2*sz;
#            diag(lind+1,4)      = -2*sx;
#            diag(lind+N.nx,5)   = -sz;
#        otherwise
#            error('Thermal boundary condition not defined')
#    end
#    switch B.rtbc
#        case {'const','dirichlet'}
#            diag(rind,3)        = 1;
#        case {'flux','neumann'}
#            diag(rind-N.nx,1)   = -sz;
#            diag(rind-1,2)      = -2*sx;
#            diag(rind,3)        = 1+2*sx+2*sz;
#            diag(rind+N.nx,5)   = -sz;
#        otherwise
#            error('Thermal boundary condition not defined')
#    end
#    
#    N.A         = spdiags(diag,[-N.nx,-1,0,1,N.nx],N.nx*N.nz,N.nx*N.nz);
#    
#    N.beenhere  = N.beenhere + 1;
#end
#
### Solve system of equations -------------------------------------------- #
#switch B.ttbc
#    case {'flux','neumann'}
#        T0(N.nz,2:(N.nx-1)) = T0(N.nz,2:(N.nx-1)) - 2*sz*N.dz*B.thf;
#    case {'const','dirichlet'}
#        Q(N.nz,:)           = 0;
#end
#switch B.btbc
#    case {'flux','neumann'}
#        T0(1,2:(N.nx-1))    = T0(1,2:(N.nx-1)) - 2*sz*N.dz*B.bhf;
#    case {'const','dirichlet'}
#        Q(1,:)              = 0;
#end
#switch B.ltbc
#    case {'flux','neumann'}
#        T0(2:(N.nz-1),1)    = T0(2:(N.nz-1),1) + 2*sx*N.dx*B.lhf;
#    case {'const','dirichlet'}
#        Q(2:N.nz-1,1)       = 0;
#end
#switch B.rtbc
#    case {'flux','neumann'}
#        T0(2:(N.nz-1),N.nx) = T0(2:(N.nz-1),N.nx) + 2*sx*N.dx*B.rhf;
#    case {'const','dirichlet'}
#        Q(2:N.nz-1,N.nx)    = 0;
#end
#
#rhsT    = reshape(T0',[N.nx*N.nz,1]);
#rhsQ    = reshape(Q',[N.nx*N.nz,1]).*dt/P.rho/P.cp;
#
#rhs     = rhsT + rhsQ;
#
## tic;
#T1          = N.A\rhs;
## tend        = toc; 
## tend
#
#T1      = reshape(T1,[N.nx,N.nz])';
#
#switch lower(B.btbc)
#    case {'flux','neumann'}
#        T1(1,1)         = T1(1,2);
#        T1(1,N.nx)      = T1(1,N.nx-1);
#end
#switch lower(B.ttbc)
#    case {'flux','neumann'}
#        T1(N.nz,1)       = T1(N.nz,2);
#        T1(N.nz,N.nx)    = T1(N.nz,N.nx-1);
#end
#
## end
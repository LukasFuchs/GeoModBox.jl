#function CNA!()
#    # dT/dt = kappa*d^2T_ij/dx_i^2 + Q_ij/rho/cp
#    # ----------------------------------------------------------------------- #
#    
#    ## Define Coeffizients for Matrix A and boundary conditions ------------- #
#    # Setup coefficient matrix A -------------------------------------------- #
#    sx      =   P.kappa*dt/N.dx^2/2;
#    sz      =   P.kappa*dt/N.dz^2/2;
#    
#    if N.beenhere == 0
#        # Erstellung der durchlaufenden Indizes ----------------------------- #
#        Number  =   zeros(N.nz,N.nx);
#        num     =   1;
#        for i=1:N.nz
#            for j=1:N.nx
#                Number(i,j) =   num;
#                num = num+1;
#            end
#        end
#        
#        # Define diagonals for matrix
#        diagL   =   zeros(N.nx*N.nz,5);
#        diagR   =   zeros(N.nx*N.nz,5);
#        
#        # Inner index ------------------------------------------------------- #
#        iind    =    Number(2:(N.nz-1),2:(N.nx-1));
#        
#        diagL(iind-N.nx,1)  =   -sz;
#        diagL(iind-1,2)     =   -sx;
#        diagL(iind,3)       =   1+2*sx+2*sz;
#        diagL(iind+1,4)     =   -sx;
#        diagL(iind+N.nx,5)  =   -sz;
#        
#        diagR(iind-N.nx,1)  =   sz;
#        diagR(iind-1,2)     =   sx;
#        diagR(iind,3)       =   1-2*sx-2*sz;
#        diagR(iind+1,4)     =   sx;
#        diagR(iind+N.nx,5)  =   sz;
#        
#        # Top and bottom index and boundary conditions ---------------------- #    
#        switch B.ttbc
#            case {'const','dirichlet'}            
#                tind                =   Number(N.nz,:);
#                diagL(tind,3)       =   1;
#                diagR(tind,3)       =   1;
#            case {'flux','neumann'}    
#                tind                =   Number(N.nz,2:(N.nx-1));
#                diagL(tind-N.nx,1)  =   -2*sz;
#                diagL(tind-1,2)     =   -sx;
#                diagL(tind,3)       =   1+2*sx+2*sz;
#                diagL(tind+1,4)     =   -sx;
#                
#                diagR(tind-N.nx,1)  =   2*sz;
#                diagR(tind-1,2)     =   sx;
#                diagR(tind,3)       =   1-2*sx-2*sz;
#                diagR(tind+1,4)     =   sx;
#                
#                diagL(N.nx*(N.nz-1)+1,3) =   1; # Left Corner grid point
#                diagL(N.nx*N.nz,3)       =   1; # Right Corner grid point
#                diagR(N.nx*(N.nz-1)+1,3) =   1; # Left Corner grid point
#                diagR(N.nx*N.nz,3)       =   1; # Right Corner grid point
#            otherwise
#                error('Thermal boundary condition not defined!')
#        end    
#        switch B.btbc
#            case {'const','dirichlet'}      
#                bind                =   Number(1,:);
#                diagL(bind,3)       =   1;
#                diagR(bind,3)       =   1;
#            case {'flux','neumann'}            
#                bind                =   Number(1,2:(N.nx-1));
#                diagL(bind-1,2)     =   -sx;
#                diagL(bind,3)       =   1+2*sx+2*sz;
#                diagL(bind+1,4)     =   -sx;
#                diagL(bind+N.nx,5)  =   -2*sz;
#                
#                diagR(bind-1,2)     =   sx;
#                diagR(bind,3)       =   1-2*sx-2*sz;
#                diagR(bind+1,4)     =   sx;
#                diagR(bind+N.nx,5)  =   2*sz;
#                
#                diagL(1,3)          =   1; # Left Corner grid point
#                diagL(N.nx,3)       =   1; # Right Corner grid point
#                diagR(1,3)          =   1; # Left Corner grid point
#                diagR(N.nx,3)       =   1; # Right Corner grid point
#            otherwise
#                error('Thermal boundary condition not defined!')
#        end
#        
#        # Left and right index and boundary conditions ---------------------- #    
#        lind                = Number(2:(N.nz-1),1);
#        rind                = Number(2:(N.nz-1),N.nx);
#        
#        switch B.ltbc
#            case {'const','dirichlet'}            
#                diagL(lind,3)       = 1;
#                diagR(lind,3)       = 1;
#            case {'flux','neumann'}            
#                diagL(lind-N.nx,1)  = -sz;
#                diagL(lind,3)       = 1+2*sx+2*sz;
#                diagL(lind+1,4)     = -2*sx;
#                diagL(lind+N.nx,5)  = -sz;
#                
#                diagR(lind-N.nx,1)  = sz;
#                diagR(lind,3)       = 1-2*sx-2*sz;
#                diagR(lind+1,4)     = 2*sx;
#                diagR(lind+N.nx,5)  = sz;
#            otherwise
#                error('Thermal boundary condition not defined')
#        end
#        switch B.rtbc
#            case {'const','dirichlet'}            
#                diagL(rind,3)       = 1;
#                diagR(rind,3)       = 1;
#            case {'flux','neumann'}
#                diagL(rind-N.nx,1)  = -sz;
#                diagL(rind-1,2)     = -2*sx;
#                diagL(rind,3)       = 1+2*sx+2*sz;
#                diagL(rind+N.nx,5)  = -sz;
#                
#                diagR(rind-N.nx,1)  = sz;
#                diagR(rind-1,2)     = 2*sx;
#                diagR(rind,3)       = 1-2*sx-2*sz;
#                diagR(rind+N.nx,5)  = sz;
#            otherwise
#                error('Thermal boundary condition not defined')
#        end
#        
#        N.AL        = spdiags(diagL,[-N.nx,-1,0,1,N.nx],N.nx*N.nz,N.nx*N.nz);
#        N.AR        = spdiags(diagR,[-N.nx,-1,0,1,N.nx],N.nx*N.nz,N.nx*N.nz);
#        
#        N.beenhere  = N.beenhere + 1;
#    end
#    
#    ## Solve system of equations -------------------------------------------- #
#    switch B.ttbc
#        case {'flux','neumann'}
#            T0(N.nz,2:(N.nx-1)) = T0(N.nz,2:(N.nx-1)) + 4*sz*N.dz*B.thf + ...
#                Q(N.nz,2:(N.nx-1)).*dt/P.rho/P.cp;
#        case {'const','dirichlet'}
#            Q(N.nz,:)           = 0;
#    end
#    switch B.btbc
#        case {'flux','neumann'}
#            T0(1,2:(N.nx-1))    = T0(1,2:(N.nx-1)) - 4*sz*N.dz*B.bhf + ...
#                Q(1,2:(N.nx-1)).*dt/P.rho/P.cp;
#        case {'const','dirichlet'}
#            Q(1,:)              = 0;
#    end
#    switch B.ltbc
#        case {'flux','neumann'}
#            T0(2:(N.nz-1),1)    = T0(2:(N.nz-1),1) - 4*sx*N.dx*B.lhf + ...
#                Q(2:N.nz-1,1).*dt/P.rho/P.cp; 
#        case {'const','dirichlet'}
#            Q(2:N.nz-1,1)       = 0;
#    end
#    switch B.rtbc
#        case {'flux','neumann'}
#            T0(2:(N.nz-1),N.nx) = T0(2:(N.nz-1),N.nx) + 4*sx*N.dx*B.rhf + ...
#                Q(2:N.nz-1,N.nx).*dt/P.rho/P.cp;
#        case {'const','dirichlet'}
#            Q(2:N.nz-1,N.nx)    = 0;
#    end
#    
#    rhsT    = reshape(T0',[N.nx*N.nz,1]);
#    rhsQ    = reshape(Q',[N.nx*N.nz,1]).*dt/P.rho/P.cp;
#    
#    rhs     = N.AR*rhsT + rhsQ;
#    
#    T1      = N.AL\rhs;
#    
#    T1      = reshape(T1,[N.nx,N.nz])';
#    
#    switch lower(B.btbc)
#        case {'flux','neumann'}
#            T1(1,1)         = T1(1,2);
#            T1(1,N.nx)      = T1(1,N.nx-1);
#    end
#    switch lower(B.ttbc)
#        case {'flux','neumann'}
#            T1(N.nz,1)       = T1(N.nz,2);
#            T1(N.nz,N.nx)    = T1(N.nz,N.nx-1);
#    end
#    
#    end
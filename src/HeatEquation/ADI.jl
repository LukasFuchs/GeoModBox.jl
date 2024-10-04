#function ADI!()
#    # Function to solve 2D heat diffusion equation using the alternating direct
#    # implicit finite difference scheme.
#    # assuming constant k, rho, cp
#    # dT/dt = kappa*d^2T_ij/dx_i^2 + Q_ij/rho/cp
#    # ----------------------------------------------------------------------- #
#    
#    ## Define Coeffizients for Matrix A and boundary conditions ------------- #
#    
#    if N.beenhere == 0
#        # Erstellung der durchlaufenden Indizes ----------------------------- #
#        # Gleichungssystem fuer ADI Solver:
#        # Erster ADI Schritt: A*T^(l+1/2) = B*T^l
#        # Zweiter ADI Schritt: C*T^(l+1) = D*T^(l+1/2)
#        Number1 = zeros(N.nz,N.nx);     # Index fuer ersten ADI Durchlauf
#        Number2 = zeros(N.nz,N.nx);     # Index fuer zweiten ADI Durchlauf
#        num     = 1;
#        for j=1:N.nx
#            for i=1:N.nz
#                # Vertikale Laufrichtung
#                Number1(i,j) = num;
#                num = num+1;
#            end
#        end
#        num     = 1;
#        for i=1:N.nz
#            for j=1:N.nx
#                # Horizontale Laufrichtung
#                Number2(i,j) = num;
#                num = num+1;
#            end
#        end
#        
#        # Setup coefficient matrices ---------------------------------------- #
#        a       = P.kappa*dt/N.dz^2/2;
#        b       = P.kappa*dt/N.dx^2/2;
#        
#        # Define diagonals for matrix
#        diag1   = zeros(N.nx*N.nz,3);
#        diag2   = zeros(N.nx*N.nz,3);
#        diag3   = zeros(N.nx*N.nz,3);
#        diag4   = zeros(N.nx*N.nz,3);
#        diag5   = zeros(N.nx*N.nz,1);
#        diag6   = zeros(N.nx*N.nz,1);
#        
#        # Inner index ------------------------------------------------------- #
#        iind1   = Number1(2:(N.nz-1),2:(N.nx-1));
#        iind2   = Number2(2:(N.nz-1),2:(N.nx-1));
#        
#        # Diagonalen fuer Koeffizientenmatrix A auf der linken Seite
#        # Vertikale Laufrichtung
#        diag1(iind1-1,1)    = -a;       # j,i-1
#        diag1(iind1,2)      = 1+2*a;    # j,i
#        diag1(iind1+1,3)    = -a;       # j,i+1
#        
#        # Diagonalen fuer Koeffizientenmatrix B auf der rechten Seite
#        # Vertikale Laufrichtung
#        diag2(iind1-N.nz,1) = b;        # j-1,i
#        diag2(iind1,2)      = 1-2*b;    # j,i
#        diag2(iind1+N.nz,3) = b;        # j+1,1
#        
#        # Diagonalen fuer Koeffizientenmatrix C auf der linke Seite
#        # Horizontale Laufrichtung
#        diag3(iind2-1,1)    = -b;       # j-1,i
#        diag3(iind2,2)      = 1+2*b;    # j,i
#        diag3(iind2+1,3)    = -b;       # j+1,i
#        
#        # Diagonalen fuer Koeffizientenmatrix D auf der rechten Seite
#        # Horizontale Laufrichtung
#        diag4(iind2-N.nx,1) = a;        # j,i-1
#        diag4(iind2,2)      = 1-2*a;    # j,i
#        diag4(iind2+N.nx,3) = a;        # j,i+1
#        
#        # Top and bottom index and boundary conditions ---------------------- #    
#        switch lower(B.ttbc)
#            case {'const','dirichlet'}
#                tind    = Number1(N.nz,:);
#                tind2   = Number2(N.nz,:);
#                
#                diag1(tind,2)   = 1;
#                diag2(tind,2)   = 1;
#                diag3(tind2,2)  = 1;
#                diag4(tind2,2)  = 1;
#            case {'flux','neumann'}
#                tind            = Number1(N.nz,2:N.nx-1);
#                tind2           = Number2(N.nz,2:N.nx-1);
#                
#                # Vertikale Laufrichtung - 1. ADI Schritt
#                diag1(N.nz,2)       = 1;        # Left lower Corner Grid point
#                diag1(tind-1,1)     = -2*a;     # j,nz-1
#                diag1(tind,2)       = 1+2*a;    # j,nz
#                diag1(N.nx*N.nz,2)  = 1;        # right lower Corner Grid point
#                
#                diag2(N.nz,2)       = 1;
#                diag2(tind-N.nz,1)  = b;        # j-1,nz
#                diag2(tind,2)       = 1-2*b;    # j,nz
#                diag2(tind+N.nz,3)  = b;        # j+1,nz
#                diag2(N.nx*N.nz,2)  = 1;
#                
#                diag5(tind,1)       = -2*a*N.dz*B.thf;
#                # Assume B.lhf is less then zero for temperature influx!
#                
#                # Horizontale Laufrichtung - 2. ADI Schritt
#                diag3(N.nx*(N.nz-1)+1,2)    = 1;
#                diag3(tind2-1,1)            = -b;       # j-1,nz
#                diag3(tind2,2)              = 1+2*b;    # j,nz
#                diag3(tind2+1,3)            = -b;       # j+1,nz
#                diag3(N.nx*N.nz,2)          = 1;
#                
#                diag4(N.nx*(N.nz-1)+1,2)    = 1;
#                diag4(tind2-N.nx,1)         = 2*a;      # 1,nz-1
#                diag4(tind2,2)              = 1-2*a;    # j,nz
#                diag4(N.nx*N.nz,2)          = 1;
#                
#                diag6(tind2,1)              = -2*a*N.dz*B.thf;
#                # Assume B.lhf is less then zero for temperature influx!
#                
#            otherwise
#                error('Thermal boundary condition not defined!')
#        end
#        switch lower(B.btbc)
#            # i = 1
#            case {'const','dirichlet'}
#                bind            = Number1(1,:);
#                bind2           = Number2(1,:);
#                
#                diag1(bind,2)   = 1;
#                diag2(bind,2)   = 1;
#                diag3(bind2,2)  = 1;
#                diag4(bind2,2)  = 1;
#                
#            case {'flux','neumann'}
#                bind            = Number1(1,2:N.nx-1);
#                bind2           = Number2(1,2:N.nx-1);
#                
#                # Vertikale Laufrichtung - 1. ADI Schritt
#                diag1(1,2)          = 1;        # Left lower Corner Grid point
#                diag1(bind,2)       = 1+2*a;    # j,1
#                diag1(bind+1,3)     = -2*a;     # j,2
#                diag1((N.nx-1)*N.nz+1,2)  = 1;  # right lower Corner Grid point
#                
#                diag2(1,2)          = 1;
#                diag2(bind-N.nz,1)  = b;        # j-1,1
#                diag2(bind,2)       = 1-2*b;    # j,1
#                diag2(bind+N.nz,3)  = b;        # j+1,1
#                diag2((N.nx-1)*N.nz+1,2)  = 1;
#                
#                diag5(bind,1)       = -2*a*N.dz*B.bhf;
#                # Assume B.lhf is less then zero for temperature influx!
#                
#                # Horizontale Laufrichtung - 2. ADI Schritt
#                diag3(1,2)          = 1;
#                diag3(bind2-1,1)    = -b;       # j-1,1
#                diag3(bind2,2)      = 1+2*b;    # j,1
#                diag3(bind2+1,3)    = -b;       # j+1,1
#                diag3(N.nx,2)       = 1;
#                
#                diag4(1,2)          = 1;
#                diag4(bind2,2)      = 1-2*a;    # j,1
#                diag4(bind2+N.nx,3) = 2*a;      # 1,i+1
#                diag4(N.nx,2)       = 1;
#                
#                diag6(bind2,1)      = -2*a*N.dz*B.bhf;
#                # Assume B.lhf is less then zero for temperature influx!
#                
#            otherwise
#                error('Thermal boundary condition not defined!')
#        end
#        
#        # Left and right index and boundary conditions ---------------------- #
#        lind    = Number1(2:(N.nz-1),1);
#        rind    = Number1(2:(N.nz-1),N.nx);
#        lind2   = Number2(2:(N.nz-1),1);
#        rind2   = Number2(2:(N.nz-1),N.nx);
#        
#        switch lower(B.ltbc)
#            case {'const','dirichlet'}
#                diag1(lind,2)    = 1;
#                diag2(lind,2)    = 1;
#                diag3(lind2,2)   = 1;
#                diag4(lind2,2)   = 1;
#                
#            case {'flux','neumann'}
#                # Vertikale Laufrichtung - 1. ADI Schritt
#                diag1(lind-1,1)     = -a;       # 1,i-1
#                diag1(lind,2)       = 1+2*a;    # 1,i
#                diag1(lind+1,3)     = -a;       # 1,i+1
#                
#                diag2(lind,2)       = 1-2*b;    # 1,i
#                diag2(lind+N.nz,3)  = 2*b;      # 2,i
#                
#                diag5(lind,1)       = 2*b*N.dx*B.lhf;
#                
#                # Horizontale Laufrichtung - 2. ADI Schritt
#                diag3(lind2,2)      = 1+2*b;    # 1,i
#                diag3(lind2+1,3)    = -2*b;     # 2,i
#                
#                diag4(lind2-N.nx,1) = a;        # 1,i-1
#                diag4(lind2,2)      = 1-2*a;    # 1,i
#                diag4(lind2+N.nx,3) = a;        # 1,i+1
#                
#                diag6(lind2,1)      = 2*b*N.dx*B.lhf;
#                
#            otherwise
#                error('Thermal boundary condition not defined')
#        end
#        switch lower(B.rtbc)
#            case {'const','dirichlet'}
#                diag1(rind,2)    = 1;
#                diag2(rind,2)    = 1;
#                diag3(rind2,2)   = 1;
#                diag4(rind2,2)   = 1;
#                
#            case {'flux','neumann'}
#                # Vertikale Laufrichtung - 1. ADI Schritt
#                diag1(rind-1,1)     = -a;       # i-1,nx
#                diag1(rind,2)       = 1+2*a;    # i,nx
#                diag1(rind+1,3)     = -a;       # i+1,nx
#                
#                diag2(rind-N.nz,1)  = 2*b;      # i,nx-1
#                diag2(rind,2)       = 1-2*b;    # i,nx
#                
#                diag5(rind,1)       = 2*b*N.dx*B.rhf;
#                
#                # Horizontale Laufrichtung - 2. ADI Schritt
#                diag3(rind2-1,1)    = -2*b;     # i,nx
#                diag3(rind2,2)      = 1+2*b;    # i,nx-1
#                
#                diag4(rind2-N.nx,1) = a;        # i-1,nx
#                diag4(rind2,2)      = 1-2*a;    # i,nx
#                diag4(rind2+N.nx,3) = a;        # i+1,nx
#                
#                diag6(rind2,1)      = 2*b*N.dx*B.rhf;
#                
#            otherwise
#                error('Thermal boundary condition not defined')
#        end
#        
#        N.A1        = spdiags(diag1,[-1,0,1],N.nx*N.nz,N.nx*N.nz);
#        N.A2        = spdiags(diag2,[-N.nz,0,N.nz],N.nx*N.nz,N.nx*N.nz);
#        N.A3        = spdiags(diag3,[-1,0,1],N.nx*N.nz,N.nx*N.nz);
#        N.A4        = spdiags(diag4,[-N.nx,0,N.nx],N.nx*N.nz,N.nx*N.nz);
#        N.Fl1       = diag5;
#        N.Fl2       = diag6;
#        
#        N.beenhere  = N.beenhere + 1;
#    end
#    
#    ## Solve system of equations -------------------------------------------- #
#    switch lower(B.ttbc)
#        case {'const','dirichlet'}
#            Q(N.nz,:)           = 0;
#    end
#    switch lower(B.btbc)
#        case {'const','dirichlet'}
#            Q(1,:)              = 0;
#    end
#    switch lower(B.ltbc)
#        case {'const','dirichlet'}
#            Q(2:N.nz-1,1)       = 0;
#    end
#    switch lower(B.rtbc)
#        case {'const','dirichlet'}
#            Q(2:N.nz-1,N.nx)    = 0;
#    end
#    
#    
#    # First ADI step ========== Vertikale Laufrichtung
#    rhsT1   = reshape(T0,[N.nx*N.nz,1]);
#    rhsQ1   = reshape(Q,[N.nx*N.nz,1]).*dt/2/P.rho/P.cp;
#    rhsT1   = N.A2*rhsT1 + rhsQ1 + N.Fl1;
#    
#    # T(l+1/2)
#    T1      = N.A1\rhsT1;
#    T1      = reshape(T1,[N.nz,N.nx]);
#    switch lower(B.btbc)
#        case {'flux','neumann'}
#            T1(1,1)         = T1(1,2);
#            T1(1,N.nx)      = T1(1,N.nx-1);
#    end
#    switch lower(B.ttbc)
#        case {'flux','neumann'}
#            T1(N.nz,1)      = T1(N.nz,2);
#            T1(N.nz,N.nx)   = T1(N.nz,N.nx-1);
#    end
#    
#    # Second ADI step ========== Horizontala Laufrichtung
#    rhsT2   = reshape(T1',[N.nx*N.nz,1]);
#    rhsQ2   = reshape(Q',[N.nx*N.nz,1]).*dt/2/P.rho/P.cp;
#    rhsT2   = N.A4*rhsT2 + rhsQ2 + N.Fl2;
#    
#    # T(l+1)
#    T2      = N.A3\rhsT2;
#    T2      = reshape(T2,[N.nx,N.nz])';
#    
#    switch lower(B.btbc)
#        case {'flux','neumann'}
#            T2(1,1)         = T2(1,2);
#            T2(1,N.nx)      = T2(1,N.nx-1);
#    end
#    switch lower(B.ttbc)
#        case {'flux','neumann'}
#            T2(N.nz,1)      = T2(N.nz,2);
#            T2(N.nz,N.nx)   = T2(N.nz,N.nx-1);
#    end
#    
#    
#    end    
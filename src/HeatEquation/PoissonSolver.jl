function [F] = SolvePoisson2D(Q,P,N,B)
    % Function to solve 2D heat diffusion equation using the explicit finite
    % difference scheme
    % [Q] = W/m^3
    % ----------------------------------------------------------------------- %
    
    %% Define Coeffizients for Matrix A and boundary conditions ------------- %
    
    
    % Erstellung der durchlaufenden Indizes ----------------------------- %
    Number  = zeros(N.nz,N.nx);
    Number2 = zeros(N.nz,N.nx);
    num     = 1;
    for i=1:N.nz
        for j=1:N.nx
            l   = j+(i-1)*N.nx;
            Number(i,j) = num;
            Number2(i,j) = l;
            num = num+1;
        end
    end
    
    % Setup coefficient matrix A ---------------------------------------- %
    a       = 1/N.dz^2;
    b       = 1/N.dx^2;
    c       = -2*(a+b);
    
    % Define diagonals for matrix
    diag    = zeros(N.nx*N.nz,5);
    
    % Inner index ------------------------------------------------------- %
    iind   = Number(2:(N.nz-1),2:(N.nx-1));
    
    diag(iind-N.nx,1)   = a;
    diag(iind-1,2)      = b;
    diag(iind,3)        = c;
    diag(iind+1,4)      = b;
    diag(iind+N.nx,5)   = a;
    
    % Top and bottom index and boundary conditions ---------------------- %
    switch B.btbc
        case {'const','dirichlet'}
            bind                = Number(1,:);
            diag(bind,3)        = 1;
        case {'flux','neumann'}
            bind                = Number(1,2:N.nx-1);        
            diag(bind-1,2)      = b;
            diag(bind,3)        = c;
            diag(bind+1,4)      = b;        
            diag(bind+N.nx,5)   = 2*a;   
            
            % Left corner
            blcind              = Number(1,1);        
            diag(blcind,3)      = c;
            diag(blcind+1,4)    = 2*b;        
            diag(blcind+N.nx,5) = 2*a; 
            
            % Right corner
            brcind              = Number(1,N.nx);         
            diag(brcind-1,2)    = 2*b; 
            diag(brcind,3)      = c;   
            diag(blcind+N.nx,5) = 2*a; 
        otherwise
            error('Thermal boundary condition not defined!')
    end
    switch B.ttbc
        case {'const','dirichlet'}
            tind                = Number(N.nz,:);
            diag(tind,3)        = 1;
        case {'flux','neumann'}
            tind                = Number(N.nz,2:N.nx-1);
            diag(tind-N.nx,1)   = 2*a;
            diag(tind-1,2)      = b;
            diag(tind,3)        = c;
            diag(tind+1,4)      = b;        
            
            % Left corner
            tlcind              = Number(N.nz,1); 
            diag(tlcind-N.nx,1) = 2*a;
            diag(tlcind,3)      = c;
            diag(tlcind+1,4)    = 2*b;
            
            % Right corner
            trcind              = Number(N.nz,N.nx); 
            diag(trcind-N.nx,1) = 2*a;
            diag(trcind-1,2)    = 2*b;
            diag(trcind,3)      = c;        
        otherwise
            error('Thermal boundary condition not defined!')
    end
    
    % Left and right index and boundary conditions ---------------------- %
    lind    = Number(2:(N.nz-1),1);
    rind    = Number(2:(N.nz-1),N.nx);
    
    switch B.ltbc
        case {'const','dirichlet'}
            diag(lind,3)    = 1;
        case {'flux','neumann'}
            diag(lind-N.nx,1)   = a;
            diag(lind,3)        = c;
            diag(lind+1,4)      = 2*b;
            diag(lind+N.nx,5)   = a;
        otherwise
            error('Thermal boundary condition not defined')
    end
    switch B.rtbc
        case {'const','dirichlet'}
            diag(rind,3)    = 1;
        case {'flux','neumann'}
            diag(rind-N.nx,1)   = a;
            diag(rind-1,2)      = 2*b;
            diag(rind,3)        = c;
            diag(rind+N.nx,5)   = a;        
        otherwise
            error('Thermal boundary condition not defined')
    end
    
    A       = spdiags(diag,[-N.nx,-1,0,1,N.nx],N.nx*N.nz,N.nx*N.nz);
    
    %% Solve system of equations -------------------------------------------- %
    rhs    = -1.*reshape(Q',[N.nx*N.nz,1])./P.k;
    
    switch B.btbc
        case {'flux','neumann'}
            %i=1
            rhs(bind)   = rhs(bind) - 2*N.dz*B.bhf*a;        
            rhs(blcind) = rhs(blcind) - 2*N.dx*B.lhf*b - 2*N.dz*B.bhf*a; 
            rhs(brcind) = rhs(brcind) - 2*N.dx*B.rhf*b + 2*N.dz*B.bhf*a; 
    end
    switch B.ttbc
        case {'flux','neumann'}
            %i=nz
            rhs(tind)   = rhs(tind) - 2*N.dz*B.thf*a;
            rhs(tlcind) = rhs(tlcind) - 2*N.dx*B.lhf*b - 2*N.dz*B.thf*a; 
            rhs(trcind) = rhs(trcind) - 2*N.dx*B.rhf*b - 2*N.dz*B.thf*a; 
    end
    switch B.ltbc
        case {'flux','neumann'}
            rhs(lind)   = rhs(lind) - 2*N.dx*B.lhf*b;
    end
    switch B.rtbc
        case {'flux','neumann'}
            rhs(rind)   = rhs(rind) - 2*N.dx*B.rhf*b;
    end
    
    F      = A\rhs;
    
    F      = reshape(F,[N.nx,N.nz])';
    
    end
#function ForwarEuler!()
#    # Function to solve 2D heat diffusion equation using the explicit finite
#    # difference scheme
#    # Q - Waermeproduktionsrate pro Volumen [W/m^3]
#    # ----------------------------------------------------------------------- #
#    
#    T1      = zeros(N.nz,N.nx);
#    
#    ind1    = 2:N.nz-1;
#    ind2    = 2:N.nx-1;
#    
#    sx      = P.kappa*dt/N.dx^2;
#    sz      = P.kappa*dt/N.dz^2;
#    
#    T1(ind1,ind2) = T0(ind1,ind2) + ...
#        sx.*(T0(ind1,ind2+1) - 2.*T0(ind1,ind2) + T0(ind1,ind2-1)) + ...
#        sz.*(T0(ind1+1,ind2) - 2.*T0(ind1,ind2) + T0(ind1-1,ind2)) + ...
#        Q(ind1,ind2)*dt/P.rho/P.cp;
#    
#    # Boundary condition ---------------------------------------------------- #
#    # Bottom ---------------------------------------------------------------- #
#    switch B.btbc
#        case {'const','dirichlet'}
#            T1(1,:)         = T0(1,:);
#        case {'flux','neumann'}
#            T1(1,ind2) = T0(1,ind2) + ...
#                sx.*(T0(1,ind2+1) - 2.*T0(1,ind2) + T0(1,ind2-1)) + ...
#                2*sz.*(T0(2,ind2) - T0(1,ind2) + N.dz*B.bhf) + ...
#                Q(1,ind2)*dt/P.rho/P.cp;
#            T1(1,1)    = T1(1,2);
#            T1(1,N.nx) = T1(1,N.nx-1);        
#    end
#    # Top ------------------------------------------------------------------- #
#    switch B.ttbc
#        case {'const','dirichlet'}
#            T1(N.nz,:)    = T0(N.nz,:);
#        case {'flux','neumann'}
#            T1(N.nz,ind2) = T0(N.nz,ind2) + ...
#                sx.*(T0(N.nz,ind2+1) - 2.*T0(N.nz,ind2) + T0(N.nz,ind2-1)) + ...
#                2*sz.*(T0(N.nz-1,ind2) - T0(N.nz,ind2) - N.dz*B.thf) + ...
#                Q(N.nz,ind2)*dt/P.rho/P.cp;
#            T1(N.nz,1)    = T1(N.nz,2);
#            T1(N.nz,N.nx) = T1(N.nz,N.nx-1);        
#    end
#    # Left ------------------------------------------------------------------ #
#    switch B.ltbc
#        case {'const','dirichlet'}
#            T1(ind1,1)         = T0(ind1,1);
#        case {'flux','neumann'}
#            T1(ind1,1)  = T0(ind1,1)  + ...
#                2*sx.*(T0(ind1,2) - T0(ind1,1) + N.dx*B.lhf) +...
#                sz.*(T0(ind1+1,1) - 2.*T0(ind1,1) + T0((ind1)-1,1)) + ...
#                Q(ind1,1)*dt/P.rho/P.cp;        
#    end
#    # Right ----------------------------------------------------------------- #
#    switch B.rtbc
#        case {'const','dirichlet'}
#            T1(ind1,N.nx)        = T0(ind1,N.nx);
#        case {'flux','neumann'}
#            T1(ind1,N.nx) = T0(ind1,N.nx) + ...
#                2*sx.*(T0(ind1,N.nx-1) - T0(ind1,N.nx) + N.dx*B.rhf) + ...
#                sz.*(T0((ind1)+1,1) - 2.*T0(ind1,1) + T0((ind1)-1,1)) + ...
#                Q(ind1,N.nx)*dt/P.rho/P.cp;
#    end
#    end
#
clear
clc

## Physikalischer Parameter --------------------------------------------- #
P.L         =   4e3;        #   [m]
P.H         =   2e3;        #   [m]

P.k1        =   5.6;        #   Waermeleitfaehigkeit, W/m/K
P.k2        =   6.2;        #   Waermeleitfaehigkeit, W/m/K

P.Wcave     =   200;        #
P.Hcave     =   200;        #
P.Dcave     =   1e3;        # Tiefe des zentrums der welle
P.Zcave     =   P.L/2; 

P.Q0        =   0.3;        # W/m³ Q = rho*H

## Numerische Parameter ------------------------------------------------- #
N.nx        =   641;        #   # Gitterpunkte in x-Richtung
N.nz        =   321;         #   # Gitterpunkte in z-Richtung

N.dx        =   P.L/(N.nx-1);   #   Gitterabstand in x-Richtung
N.dz        =   P.H/(N.nz-1);   #   Gitterabstand in z-Richtung

B.btbc      = 'const';
B.ttbc      = 'const';
B.ltbc      = 'const';
B.rtbc      = 'const';

# ------------------
B.rhf       = 0;        
B.lhf       = 0; 
# ------------------
B.bhf       = 0;        
B.thf       = 0; 

# Erstellung des Gitters ------------------------------------------------ #
[N.x,N.z]   = meshgrid(0:N.dx:P.L,-P.H:N.dz:0);


## Erstellung des Anfangstemperaturfeld --------------------------------- #
D.Q         = 0*ones(N.nz,N.nx);

P.ind2      = N.x>P.H/2; 
D.k         = P.k1.*ones(N.nz,N.nx);
D.k(P.ind2) = P.k2; 

D.T         = zeros(N.nz,N.nx);

# Defniere die Region of additional heat source
P.ind       = ((N.x>=(P.Zcave-P.Wcave/2))&(N.x<=(P.Zcave+P.Wcave/2))&...
    (N.z>=-P.Dcave-P.Hcave/2)&(N.z<=-P.Dcave+P.Hcave/2));

D.Q(P.ind)  = P.Q0;

## Solve equation ------------------------------------------------------- #
D.T         = SolvePoisson2Dvaryk(D.Q,P,N,B,D.k);

## Plot solution -------------------------------------------------------- #
set(figure(1),'Position',[331.4,231.4,914.4,420])
clf
pcolor(N.x/1e3,N.z/1e3,D.T); shading interp; colorbar
hold on
contour(N.x/1e3,N.z/1e3,D.T,100:100:1500,'k');
xlabel('x [km]'); ylabel('z [km]'); zlabel('Temperature [^oC]')
title(['Stationary temperature field for ',B.ltbc,' lateral boundary conditions'])
caxis([0 900])
set(gca,'FontWeight','Bold')















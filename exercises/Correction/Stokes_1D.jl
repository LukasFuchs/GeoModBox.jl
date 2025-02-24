clear
clc

%% Model Geometry ======================================================= %
H       =   -400e3;             % [ m ]
vx0     =   5;                  % [ cm/a ]
vx0     =   vx0/100/31536000;   % [ m/s ]
eta     =   1e21;               % [ Pa s ]
dPdx    =   -1e2;               % [ Pa/m ]; z.B. -1e2 (Warum negativ?)

%% Numerische Parameter ================================================= %
nz      =   101; 
dz      =   H/(nz-1); 

% Koordinaten
z       =   H:-dz:0;
% Index der internen Gitterpunkte
indi    =   2:nz-1; 

%% Erstellung der Koeffizientenmatrix A 
diag            =   zeros(nz,3); 

% Bestimmung der Diagonalen 
diag(indi-1,1)  =   1/dz^2; 
diag(indi,2)    =   -2/dz^2; 
diag(indi+1,3)  =   1/dz^2; 

% Randbedingungen - no slip, d.h. konstante Geschwindigkeit
diag(1,2)       =   1; 
diag(nz,2)      =   1; 

A               =   spdiags(diag,[-1 0 1],nz,nz);

%% Definition der rechten Seite des Gleichungssystems
rhs             =   zeros(nz,1); 

rhs(indi)       =   dPdx/eta;
rhs(1)          =   0; 
rhs(nz)         =   vx0; 

vx              =   A\rhs; 

%% Analytische Loesung ================================================== %
vx_ana          =   1/2/eta*dPdx.*(z.^2 - H.*z) - vx0.*z./H + vx0;

%% Darstellung der Daten
figure(1)
clf
plot(vx,z./1e3,'LineWidth',2)
hold on
plot(vx_ana,z./1e3,'r--','LineWidth',2)
xlabel('v_x [ m/s ]'); ylabel('z [km]'); title('Geschwindigkeitsprofil')
set(gca,'FontWeight','Bold','FontSize',15,'LineWidth',2);

%% ===================================================================== %%
% ============================== ENDE =================================== %
% ======================================================================= %












clear
clc

%% Model Geometry ======================================================= %
H       =   -400e3;             % [ m ]
vx0     =   2;                  % [ cm/a ]
vx0     =   vx0/100/31536000;   % [ m/s ]
dPdx    =   -1e2;               % [ Pa/m ]; z.B. -7e1 (Warum negativ?)

%% Numerische Parameter ================================================= %
nz      =   51;
dz      =   H/(nz-1);

% Koordinaten
z       =   H:-dz:0;
% Index der internen Gitterpunkte
indi    =   2:nz-1;

%% Definition der Vikosität ============================================= %
eta0        =   22;                 % Power der Viskosität oben
eta1        =   21;                 % Power der Viskosität unten
m           =   10^eta1/10^eta0;    % Viskositätsverhältnis
% eta         =   10^eta0.*m.^(z'./H);   % Viskositätsprofil
eta         =   10.^eta0*exp(log(m).*z'./H);

%% Analytische Loesung ================================================== %
vx_ana      =   -dPdx*H/10^eta0/log(m)/(m-1).*...
    (z.*(m.^((-z+H)./H) - m.^(-z./H)) + H.*(m.^(-z./H) - 1)) + ...
    vx0/(m-1).*(m.^((-z+H)./H) - 1);

%% Erstellung der Koeffizientenmatrix A
diag        =   zeros(nz,3);

% Bestimmung der Diagonalen
diag(indi-1,1)  =   (eta(indi)+eta(indi-1))/2/dz^2;
diag(indi,2)    =   -(eta(indi-1)+2*eta(indi)+eta(indi+1))/2/dz^2;
diag(indi+1,3)  =   (eta(indi)+eta(indi+1))/2/dz^2;

% Randbedingungen - no slip, d.h. konstante Geschwindigkeit
diag(1,2)       =   1;
diag(nz,2)      =   1;

A               =   spdiags(diag,[-1 0 1],nz,nz);

%% Definition der rechten Seite des Gleichungssystems
rhs             =   zeros(nz,1);

rhs(indi)       =   dPdx;
rhs(1)          =   0;          % unten
rhs(nz)         =   vx0;        % oben

vx              =   A\rhs;

%% Darstellung der Daten
figure(1)
clf
subplot(1,2,1)
plot(vx,z./1e3,'LineWidth',2)
hold on
plot(vx_ana,z./1e3,'r--','LineWidth',2)
xlabel('v_x [ m/s ]'); ylabel('z [km]'); title('Geschwindigkeitsprofil')
set(gca,'FontWeight','Bold','FontSize',15,'LineWidth',2);
subplot(1,2,2)
plot(eta,z./1e3,'LineWidth',2)
xlabel('\eta[ Pa s ]'); ylabel('z [km]'); title('Viskosity')
set(gca,'xscale','log','FontWeight','Bold','FontSize',15,'LineWidth',2);


%% ===================================================================== %%
% ============================== ENDE =================================== %
% ======================================================================= %













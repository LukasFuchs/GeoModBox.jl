function [Tnew] = SemiLagAdvection2D(vx,vz,vxo,vzo,X,Z,T,dt)
% Funktion zur Loesung der Advektionsgleichung
%
%               dA/dt = -vx * ( dA/dx ) - vz * ( dA/dz )
%
% in einer zweidimensionalen Umgebung mit Hilfe des Semi-Lagrangian Schema,
% und einer Zentralen-Punkt Iteration, d.h.
%
% ======================================================================= %

% mid-point iteration scheme -------------------------------------------- %
if isempty(vxo)
    % Im Falle das die Geschwindigkeit zeitlich konstant ist, wird die
    % aktuelle Geschwindigkeit auf die alte Geschwindigkeit
    % uebertragen.
    vxo     = vx;
    vzo     = vz;
end

% Mittlere Geschwindigkeit am Zentralen Punkt in der Zeit --------------- %
vxm     = 0.5.*(vxo + vx);
vzm     = 0.5.*(vzo + vz);

% Initialisierung der Geschwindigkeit fuer die Iteration ---------------- %
vxi     = vx;
vzi     = vz;

% Iteration ------------------------------------------------------------- %
for k = 1:10
    xp  = X - 0.5*dt.*vxi;
    zp  = Z - 0.5*dt.*vzi;
    
    vxi = interp2(X,Z,vxm,xp,zp,'linear');
    vzi = interp2(X,Z,vzm,xp,zp,'linear');
    
    vxi(isnan(vxi)) = vxm(isnan(vxi));
    vzi(isnan(vzi)) = vzm(isnan(vzi));
end
xp  = X - dt.*vxi;
zp  = Z - dt.*vzi;

Tnew    = interp2(X,Z,T,xp,zp,'cubic');

Tnew(isnan(Tnew)) = T(isnan(Tnew));

end
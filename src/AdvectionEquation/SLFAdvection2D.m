function [A] = SLFAdvection2D(vx,vz,Aold,A,dx,dz,dt)
% Funktion zur Loesung der Advektionsgleichung 
% 
%               dA/dt = -vx * ( dA/dx ) - vz * ( dA/dz )
% 
% in einer zweidimensionalen Umgebung mit Hilfe des staggered-leaped frog
% Schemas (SLF), d.h. 
% die Zeit- und Raum-Differentiale werden beide als zentrale finite
% Differenzen genommen. Daher benoetigen wir auch das ALTE Temperaturfeld! 
% 
% Das SLF Schema ist numerisch stabile fuer ein Courant Kriterium
% kleiner als 1, besitzt aber eine hoehere Fehlerordnung, 
% sowohl im Raum (O(dx^2)), als auch in der Zeit (O(dt^2)). 
% Allerdings tritt auch hier KEINE numerische Diffusion auf. 
% ======================================================================= %

nz      = size(vx,1);
nx      = size(vx,2);

% Bestimmen wir die Indizes der inneren Gitterpunkte:
iindx   = 2:(nx-1);
iindz   = 2:(nz-1);

A(iindz,iindx) = Aold(iindz,iindx) - ...
    vx(iindz,iindx).*dt/dx.*(A(iindz,iindx+1)-A(iindz,iindx-1)) - ...
    vz(iindz,iindx).*dt/dz.*(A(iindz+1,iindx)-A(iindz-1,iindx));
end
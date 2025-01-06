function [A] = UpwindAdvection2D(vx,vz,A,dx,dz,dt)
% Funktion zur Loesung der Advektionsgleichung 
% 
%               dA/dt = -vx * ( dA/dx ) - vz * ( dA/dz )
% 
% in einer zweidimensionalen Umgebung mit Hilfe des Upwind Schemas, d.h.
% die einseitige finite Differenz im Raum wird entlang der 
% Stroemungsrichtung genommen. 
% Das Upwind Schema ist numerisch stabile fuer ein Courant Kriterium
% kleiner als 1. Allerdings gibt es eine sehr starke numerische Diffusion. 
%
%======================================================================== %

nz      = size(vx,1);
nx      = size(vx,2);

% Bestimmen wir die Indizes der inneren Gitterpunkte:
iindx   = 2:(nx-1);
iindz   = 2:(nz-1);

% Ersetze die Fragezeichen gemaess den Angaben aus der Vorlesung! 

% Um zu bestimmen an welchen Gitterpunkten die Geschwindigkeit positiv
% bzw. negativ ist wenden wir einen kleinen Trick an. Dazu erzeugen wir
% jeweils 4 Arrays, in welchen die Element gleich 1 sind, wenn die
% Geschwindigkeit positiv/negativ ist und gleich null wenn das Gegenteil
% der Fall ist: 

% Setzt die Elemente gleich 1, falls vx positiv oder negativ ist. 
Dx1     = vx(iindz,iindx)>0;     Dx2     = vx(iindz,iindx)<0;

% Setzt die Elemente gleich 1, falls vz positv oder negativ ist. 
Dz1     = vz(iindz,iindx)>0;     Dz2     = vz(iindz,iindx)<0;

A(iindz,iindx) = A(iindz,iindx) - ...
    Dx1.*(vx(iindz,iindx).*dt/dx.*(A(iindz,iindx) - A(iindz,iindx-1))) - ...
    Dx2.*(vx(iindz,iindx).*dt/dx.*(A(iindz,iindx+1) - A(iindz,iindx))) - ...
    Dz1.*(vz(iindz,iindx).*dt/dz.*(A(iindz,iindx) - A(iindz-1,iindx))) - ...
    Dz2.*(vz(iindz,iindx).*dt/dz.*(A(iindz+1,iindx) - A(iindz,iindx)));

end
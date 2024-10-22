%% Advektionsgleichung 1-D
% Bisher haben wir die Energieerhaltungsgleichung nur ohne den Transport von 
% Materialien (der Advektion) betrachtet. Häufig haben wir aber ein Problem in 
% dem sich das Material bewegt und bestimmte Größen, wie z.B. die Temperatur, 
% Dichte, etc., transportiert werden müssen (z.B. in Mantleplumes). Generell gesehen, 
% ist die Mantelkonvektion ein Beispiel eines Systems, in welchem die Temperatur 
% sowohl durch *Diffusion* (vor allem in den Grenzschichten) als auch *Advektion* 
% (vor allem in Inneren) transportiert wird. 
% 
% Betrachten wir im Folgenden nur die *Advektion* alleine (d.h. _k,_ _κ_ = 0). 
% 
% Im 1-D Fall reduziert sich die Gleichung dadurch zur reinen Advektionsgleichung: 
% 
% $$\frac{\partial T}{\partial t}=-v_x \frac{\partial T}{\partial x}$$
% 
% Diese Gleichung lässt sich durch verschiedene Diskretisierungsverfahren numerisch 
% lösen und wir wollen verschiedene Verfahren programieren und auf ein bestimmtes 
% Advektionsproblem anwenden. 

% Funktion Avection1D_discretization
% ----------------------------------------------------------------------- %
% Funktion zur Lösung der 1-D Advektionsgleichung unter der Annahme einer 
% konstanten horizontalen Geschwindigkeit. 
% ----------------------------------------------------------------------- %
% Vers. 1.0 - 7.12.2020
% ======================================================================= %
clear
clc
%% 
% Betrachten wir dazu zwei 1-D Probleme (z.B. ein horizontales Temperaturprofil) 
% mit bestimmten Anomalien:
%% 
% * Einen Gaußchen Temperaturverlauf (ein glatter Übergang):
%% 
% 
%% 
% * Eine block-förmige Temperaturanomalie (ein sehr scharfer Übergang):
%% 
% 
% 
% Erstellen wir dazu erst einmal wieder unsere Geometrie und die benötigten 
% numerischen Parameter (d.h. Gitterauflösung, Gitter, etc.): 

% Geometrische Konstanten ----------------------------------------------- %
xmin    = 0;                            % [ m ]
xmax    = 40;                           % [ m ]

% Numerische Konstanten ------------------------------------------------- %
nx      = 41;                          % Anzahl der Gitterpunkte
dx      = xmax/(nx-1);                  % Gitterlänge

x       = xmin:dx:xmax;                 % x-Koordinate
xhalf   = xmax/4;

ind     = 2:nx-1;                       % Index der inneren Gitterpunkte

% Maximale Laufzeit des Models ------------------------------------------ %
tmax    = 20;                           % [ s ]

% Horizontale Geschwindigkeit ------------------------------------------- %
vx      = 1;                            % [ m/s ]

% Definition der Zeitschrittlänge --------------------------------------- %
dtfac   = 1;                          % Courant-Kriterium
dt      = dtfac*dx/abs(vx);
nt      = ceil(tmax/dt);                % Anzahl der Zeitschritte
%% 
% Nun können wir zum einen die FD-Methode wählen ('FTCS', 'upwind', 'downwind', 
% 'lax', 'slf', 'semi-lag' - die Gleichungen und Erläuterungen dazu sind im Detail 
% in den Folien der Vorlesung zu finden) und zum anderen das Anfangsprofil wählen 
% ('block' oder 'gaussian'):

fdmethod    = 'tracers';
Tanomaly    = 'block';
%% 
% Falls wir die Tracer Methode verwenden wollen, müssen wir noch die Anzahl 
% der Tracer pro Gitterlänge festlegen: 

% Tracer advection method ----------------------------------------------- %
nmx     = 15;    % Number of tracers per "cell"
%% 
% Wollen wir nun die Anfangsbedingung (d.h. das Anfangstemperaturprofil) definieren: 

switch Tanomaly
    case 'block'
        % Hintergrundtemperatur ----------------------------------------- %
        Tb      = 1000;                 % [ K ]
        
        % Lokalität und Intensität der Temperaturanomalie --------------- %
        xTl     = xmax/10;
        xTr     = xTl + xmax/10;
        Ta      = 1500;                 % [ K ]
        
        % Erstellung des Anfangstemperaturprofiles ---------------------- %
        T     = Tb.*ones(1,nx);
        T(x>=xTl&x<=xTr) = Ta;
        tc  = 100;
        
    case 'gaussian'
        % Gaußsche Temperature Verteilung ------------------------------- %
        Tb      = 1000;                 % Hintergrundtempertur
        Ampl    = 500;                  % Amplitude
        sigma   = 1;                    % Standard Abweichung
        xc      = 10;                   % x-Koordinate des Maximums
        T       = Tb + Ampl.*exp(-((x - xc).^2)./sigma^2);
        
        Tb      = min(T);
        Ta      = max(T);
        tc      = Ampl;
end

% Erstellung des Temperaturfeldes für den Ort (Spalten) und die Zeit
% (Zeilen)
T       = T.*ones(nt,nx);
%% 
% Für die Tracermethode müssen noch ein paar zusätzliche Parameter definiert 
% werden: 

switch fdmethod
    case 'tracers'
        % Gesamtanzahl der Tracer
        nm          = (nx-1)*nmx;
        % Abstand der Tracer
        dmx         = (abs(xmin)+abs(xmax))/(nm-1);
        % x-Koordinaten der Tracer
        xm          = linspace(xmin,xmax-dmx,nm) + rand(1,nm)*0.5*dmx;
        % Temperatur auf den Tracern
        Tm          = zeros(1,nm);
end
%% 
% Definieren wir nur wieder die Eigenschaften der Figur und den Namen der Datei 
% zum Speichern der Zeitentwicklung: 

% Animation settings
filename    = ['1D_Advection',num2str(dtfac),fdmethod,'.gif'];
h           = figure(1);
txt = ['Numerical solution \Deltat = ',num2str(dtfac),'\Deltax/v_x'];
%% 
% Nun können wir die Gleichungen zur Lösung der Advektionsgleichung programieren. 
% Hierbei muss man nur darauf achten, dass die Temperatur für jeden Zeitschritt 
% in dem Temperaturarray T(nt,nx) gespeichert wird, d.h. die Temperatur des aktuellen 
% Zeitpunktes befindet sich in der Zeile T(t,:).
% 
% Da wir zur Berechnung der Temperatur des neuen Zeitpunktes auf die Temperatur 
% des vorangegangenen Zeitpunktes angewiesen sind beginnen wir die Zeitschleife 
% bei t = 2, d.h. die uns bekannte Temperatur befindes sich im vorangegangenen 
% Zeitpunkt T(t-1,:). Da wir nur die Punkte im Inneren betrachten, können wir 
% das oben erstellte Index Array für die inneren Punkte (ind) verwenden, d.h die 
% Temperatur auf den inneren Punkten zum neuen Zeitpunkt befindet sich auf T(t,ind) 
% und zum uns bekannten, dem alten, Zeitpunkt auf T(t-1,ind). 

% Lösen der Advektionsgleichung --------------------------------------- %
for t=2:nt
    disp([' Time step: ',num2str(t)])
    
    switch fdmethod
        case 'FTCS'
            T(t,ind) = ...
                T(t-1,ind) - (vx*dt/2/dx)*(T(t-1,ind+1)-T(t-1,ind-1));
        case 'upwind'
            if(vx > 0)
                T(t,ind) = ...
                    T(t-1,ind) - vx*dt/dx*(T(t-1,ind)-T(t-1,ind-1));
            elseif(vx < 0)
                T(t,ind) = ...
                    T(t-1,ind) - vx*dt/dx*(T(t-1,ind+1)-T(t-1,ind));
            end
        case 'downwind'
            T(t,ind) = ...
                T(t-1,ind) - vx*dt/dx*(T(t-1,ind+1)-T(t-1,ind));
        case 'lax'
            T(t,ind) = (T(t-1,ind+1) + T(t-1,ind-1))./2 ...
                - (vx*dt/2/dx)*(T(t-1,ind+1)-T(t-1,ind-1));
        case 'slf'
            % Hierzu nehmen wir an, dass beim ersten Zeitschritt gilt:
            % T(n-1,:) = T(n,:)
            if t==2
                T(t,ind) = ...
                    T(t-1,ind) - vx*dt/dx*(T(t-1,ind+1)-T(t-1,ind-1));
            else
                T(t,ind) = ...
                    T(t-2,ind) - vx*dt/dx*(T(t-1,ind+1)-T(t-1,ind-1));
            end
        case 'semi-lag'
            % Ausgangsposition des Partikels            
            X       = x - dt.*vx;
            
            % Interpolation der Temperatur vom Gitter x auf die
            % Ausgangsposition der Partikel X
%             Tsl     = interp1(x,T(t-1,:),X,'spline','extrap');
            Tsl     = interp1(x,T(t-1,:),X,'makima');
            
            % Zuweisung der Temperatur der "advektierten" Partikel auf das 
            % reguläre Gitter: T(n+1,xi) = T(n,X)
            T(t,:)  = Tsl;
        case 'tracers'
            Tm(1,:) = interp1(x,T(t-1,:),xm,'pchip');
            xm      = AdvectMarker1D(xm,dt,vx,max(x));
            T(t,:)  = interp1(xm,Tm(1,:),x,'pchip');
    end
    
    % Mirror boundary conditions ---------------------------------------- %
    if vx > 0
        T(t,end)    = T(t,end-1); 
        T(t,1)      = T(t,end); 
    elseif vx < 0
        T(t,1)      = T(t,2); 
        T(t,end)    = T(t,1);
    end
    
    % Darstellung des Profils
    if mod(t,0)
        figure(1),clf
        plot(x,T(t,:),'o-')
        text(xhalf,Ta+tc/2,txt)
        switch fdmethod
            case 'tracers'
                hold on
                plot(xm,Tm,'.')
        end
        xlabel('x [m]'); ylabel('T [ K ]')
        title(['1-D numerical advection, scheme: ',fdmethod])
        axis([xmin xmax Tb-tc Ta+tc])
        drawnow
        
%         % Capture the plot as an image
%         frame       = getframe(h);
%         im          = frame2im(frame);
%         [imind,cm]  = rgb2ind(im,256);
%         
%         % Write to the GIF File
%         if t == 2
%             imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%         else
%             imwrite(imind,cm,filename,'gif','WriteMode','append');
%         end
    end
end
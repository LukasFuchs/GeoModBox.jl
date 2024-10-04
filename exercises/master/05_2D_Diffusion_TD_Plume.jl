## Diffusionsgleichung (2D)

clear
clc
# 'Plume', 'Anomaly', 'Gaussian', Quelle'
Model       =   'Anomaly';
# 'explicit','implicit','CNV','ADI','implicito'
FDSchema    =   'CNV';

# Physikalischer Parameter ---------------------------------------------- #
P.L         =   200e3;          #   L??nge des Models    [m]
P.H         =   200e3;          #   H??he des Models     [m]

P.k         =   6;              #   Thermische Konduktivit??t [W/m/K]
P.cp        =   1000;           #   W??rmekapazit??t [J/kg/K]
P.rho       =   3200;           #   Dichte [kg/m^3]
P.kappa     =   P.k/P.rho/P.cp; #   Thermische Diffusivit??t [m^2/s]
P.K0        =   273.15;         #   Kelvin bei 0 C

P.Q0        =   0;              #   Hintergrund Waermeproduktionsrate;

##

switch Model
    case 'Anomaly'
        T.tmax      =   200;            # Maximale Laufzeit des Models in Ma
        N.dn        =   25;             # Inkremente der Graphischen Darstellung, d.h. hier nur jeder 25 Zeitschritt
        
        P.Ttop      =   0;              # Temperatur an der Oberfl??che [ C ]
        P.Ttop      =   P.Ttop+P.K0;    # Temperatur in Kelvin
        P.Tbot      =   1300;           # Temperatur am unteren Rand  [ C ]
        P.Tbot      =   P.Tbot+P.K0;    # Temperatur in Kelvin
        
        P.Wwave     =   P.L;            # Breite der Anomaly
        P.Hwave     =   10e3;           # H??he der Anomaly
        P.Dwave     =   P.H/2;          # Tiefe des Zentrums der Anomaly
        P.Zwave     =   P.L/2;          # x-Koordinate des Zentrums der Anomaly
        
        P.Twave     =   1600;           # Temperatur der Anomalie [ C ]
        P.Twave     =   P.Twave+P.K0;   # Temperatur in Kelvin
        P.Qwave     =   2.7e-6;         # Waermeproduktionsrate pro Volumen [W/m^3]
        
    case 'Plume'
        T.tmax      =   200;            # Maximale Laufzeit des Models in Ma
        N.dn        =   25;             # Inkremente der Graphischen Darstellung, d.h. hier nur jeder 25 Zeitschritt
        
        P.Tbot      =   1300;           # Temperatur am unteren Rand  [C]
        P.Tbot      =   P.Tbot+P.K0;    # Temperatur in Kelvin
        P.Ttop      =   0;              # Temperatur am oberen Rand   [C]
        P.Ttop      =   P.Ttop+P.K0;    # Temperatur in Kelvin
        
        P.Tplume    =   2000;           # Temperatur des Plumes       [C]
        P.Tplume    =   P.Tplume+P.K0;  # Temperatur in Kelvin
        P.Wplume    =   50e3;           # Breite des Plumes [m]
        P.Xplume    =   P.L/2;          # X-Koordinate des Zentrums des Plumes
    case 'Gaussian'        
        P.L         =   200e3;          #   L??nge des Models    [m]
        P.H         =   200e3;          #   H??he des Models     [m]

        T.tmax      =   10;             # Maximale Laufzeit des Models in Ma
        N.dn        =   5;              # Inkremente der Graphischen Darstellung, d.h. hier nur jeder 25 Zeitschritt
        
        P.Tbot      =   500;            # Temperatur am unteren Rand  [C]
        P.Tbot      =   P.Tbot+P.K0;    # Temperatur in Kelvin
        P.Ttop      =   500;            # Temperatur am oberen Rand   [C]
        P.Ttop      =   P.Ttop+P.K0;    # Temperatur in Kelvin
        
        P.Tmax      =   500;            # Temperaturamplitude [C]
        P.Tmax      =   P.Tmax+P.K0;    # Temperatur in Kelvin
        P.sigma     =   8e3;            #
        P.Xc        =   P.L/2;          # X-Koordinate des Zentrums
        P.Zc        =   P.H/2;          # z-Koordinate des Zentrums
    case 'Quelle'
        P.L         =   4e3;            #   L??nge des Models    [m]
        P.H         =   2e3  ;          #   H??he des Models     [m]
        
        P.k         =   5.6;            #   Thermische Konduktivit??t [W/m/K]
        P.cp        =   1055;           #   W??rmekapazit??t [J/kg/K]
        P.rho       =   2200;           #   Dichte [kg/m^3]
        P.kappa     =   P.k/P.rho/P.cp; #   Thermische Diffusivit??t [m^2/s]
        P.K0        =   273.15;         #   Kelvin bei 0 C

        T.tmax      =   5e-3;           # Maximale Laufzeit des Models in Ma
        N.dn        =   25;             # Inkremente der Graphischen Darstellung, d.h. hier nur jeder 25 Zeitschritt
        
        P.Tback     =   0;              # Hintergrundtemperatur
        P.Tback     =   P.Tback+P.K0;   # Temperatur in Kelvin
        
        P.Wquelle   =   200;            # Breite der Quelle
        P.Hquelle   =   200;            # H??he der Quelle
        P.Dquelle   =   P.H/2;          # Tiefe des Zentrums der Quelle
        P.Zquelle   =   P.L/2;          # x-Koordinate des Zentrums der Quelle
        
        P.Qquelle   =   0.3;            # [ W/m^3 ]
end
##

# Numerische Parameter -------------------------------------------------- #
N.nx        =   51;            # # Gitterpunkte in x-Richtung
N.nz        =   51;             # # Gitterpunkte in z-Richtung
N.dx        =   P.L/(N.nx-1);   # Gitterabstand in x-Richtung
N.dz        =   P.H/(N.nz-1);   # Gitterabstand in z-Richtung

[N.x,N.z]   =   meshgrid(0:N.dx:P.L,-P.H:N.dz:0);

N.beenhere  =   0;              # Dummy parameter, nicht ver??ndern!

# Berechnung der stabilen Zeitschrittlaenge ----------------------------- #
T.day       =   3600*24;            # # Sekunden pro Tag
T.year      =   365.25*T.day;       # # Sekunden pro Jahr
T.tmax      =   T.tmax*1e6*T.year;  # Maximale Zeit in Sekunden

T.dtfac     =   1.9;                # Multiplikationsfaktor f??r dt
T.dt        =   T.dtfac*(1/(2*P.kappa*(1/N.dx^2+1/N.dz^2)));

N.nt        =   ceil(T.tmax/T.dt);  # Anzahl der Zeitschritte
T.time      =   zeros(1,N.nt);
##
B.ttbc      =   'const';            # Oben - bisher nur const
B.btbc      =   'const';            # Unten - bisher nur const
B.ltbc      =   'flux';             # links - const od. flux
B.rtbc      =   'flux';             # rechts - const od. flux
##
# Waermefluss
# ----------
B.lhf       = 0;    B.rhf       = 0;
# ----------
B.bhf       = 0;    B.thf       = 0;
## Erstellung des Anfangstemperaturfeld
# Jetzt m??ssen wir noch die Anfangsbedingungen f??r unser beider Probleme definieren:

switch Model
    case 'Plume'
        # Hintergrundfeld f??r W??rmequellen
        D.Q             =   P.Q0.*ones(N.nz,N.nx);
        
        # Temperatur der Lithosph??re - linear zunehmenden mit der Tiefe
        D.T0            =   P.Ttop + abs(N.z./P.H)*P.Tbot;
        
        # Index f??r die Gitterpunkte am unteren Rand welche ??ber dem Plume
        # liegen
        P.ind           =   find(abs(N.x(1,:)-P.Xplume) <= P.Wplume/2);
        
        # Temperatur des Plumes
        D.T0(1,P.ind)   =   P.Tplume;
        
    case 'Anomaly'
        # Hintergrundfeld f??r W??rmequellen
        D.Q             =   P.Q0.*ones(N.nz,N.nx);
        
        # Temperatur der Lithosph??re - linear zunehmenden mit der Tiefe
        D.T0            =   P.Ttop + abs(N.z./P.H)*P.Tbot;
        
        # Defniere die Region der Anomaly
        P.ind   = ((N.x>=P.Zwave-P.Wwave/2)&(N.x<=P.Zwave+P.Wwave/2)&...
            (N.z>=-P.Dwave-P.Hwave/2)&(N.z<=-P.Dwave+P.Hwave/2));
        
        # Temperatur der Anomaly
        D.T0(P.ind)     =   P.Twave;
        
        # W??rmeproduktionsrate der Anomaly
        D.Q(P.ind)      =   P.Qwave;
        
    case 'Gaussian'
        # Temperaturfeld
        D.T0    =  P.Tbot + P.Tmax .*...
                    exp( -( (N.x-P.Xc).^2 + (N.z+P.Zc).^2 ) ./...
                    ( 2*P.sigma^2/pi ));
        D.Tana          = D.T0; 
        
        # Hintergrundfeld f??r W??rmequellen
        D.Q             = P.Q0.*ones(N.nz,N.nx);
        
        D.Tmaxa         = zeros(1,N.nt); 
        D.TProfilea     = zeros(N.nz,N.nt);
        D.RMS           = zeros(1,N.nt); 
        
    case 'Quelle'
        # Hintergrundfeld f??r W??rmequellen
        D.Q             =   P.Q0.*ones(N.nz,N.nx);
        
        # Temperaturfeld
        D.T0            = ones(N.nz,N.nx).*P.Tback;
        
        # Defniere die Region der Quelle
        P.ind   = ((N.x>=P.Zquelle-P.Wquelle/2)&(N.x<=P.Zquelle+P.Wquelle/2)&...
            (N.z>=-P.Dquelle-P.Hquelle/2)&(N.z<=-P.Dquelle+P.Hquelle/2));
        
        # Waermequellen
        D.Q(P.ind)      =   P.Qquelle;
end

# Definition des Names der Abbildung die als GIF gespeichert wird.
# Animation settings
filename    = ['2D_',FDSchema,'_',Model,'.gif'];
h           = figure(1);
## Zeitschleife

D.Tmax      = zeros(1,N.nt);
D.TProfile  = zeros(N.nz,N.nt);
for n=1:N.nt
    if n>1
        switch FDSchema
            case 'explicit'
                D.T1        =   SolveDiff2Dexplicit(D.T0,D.Q,T.dt,P,N,B);                
            case 'implicit'
                [D.T1,N]    =   SolveDiff2Dimplicit(D.T0,D.Q,T.dt,P,N,B);
            case 'implicito'
                [D.T1,N]    =   SolveDiff2Dimplicit_opt(D.T0,D.Q,T.dt,P,N,B);                
            case 'ADI'
                if(n==1)
                    N.beenhere = 0;
                end
                [D.T1,N]    =   SolveDiff2DADI(D.T0,D.Q,T.dt,P,N,B);
            case 'CNV'
                if(n==1)
                    N.beenhere = 0;
                end
                [D.T1,N]    = SolveDiff2DCNV(D.T0,D.Q,T.dt,P,N,B);                
        end                            
        
        # Zuweisung der neuen Temperatur
        D.T0        =   D.T1;
        
        T.time(n)   =   T.time(n-1) + T.dt;
        
        switch Model
            case 'Gaussian'
                D.Tana      =  P.Tbot + P.Tmax ./...
                    (1 + 2*pi*T.time(n)*P.kappa/P.sigma^2) .*...
                    exp( -( (N.x-P.Xc).^2 + (N.z+P.Zc).^2 ) ./...
                    ( 2*P.sigma^2/pi + 4*T.time(n)*P.kappa));
        end    
    end
    
    switch Model
        case 'Gaussian'
            D.Tmax(n)           = max(max(D.T0)); 
            D.Tmaxa(n)          = max(max(D.Tana));
            D.TProfile(:,n)     = D.T0(N.x == P.L/2);
            D.TProfilea(:,n)    = D.Tana(N.x == P.L/2);
            
            D.epsT      = (D.Tana-D.T0);
            
            D.RMS(n)    = sqrt(sum(sum((D.epsT).^2))/(N.nx*N.nz));
        otherwise
            # Speicher das Temperaturprofil bei x = L/2
            D.TProfile(:,n)     =   D.T0(N.x == P.L/2);
            D.Tmax(n)           =   max(D.T0(N.z == -P.H/2));
    end
    
    if (~mod(n,N.dn)||n==1||n==N.nt)
        figure(1),
        clf
        pcolor(N.x/1e3,N.z/1e3,D.T0-P.K0); shading interp; colorbar
        hold on
        contour(N.x/1e3,N.z/1e3,D.T0-P.K0,10,'k')
        switch Model
            case 'Gaussian'
                contour(N.x/1e3,N.z/1e3,D.Tana-P.K0,10,'y--')
        end
        xlabel('x [km]'); ylabel('z [km]'); zlabel('Temperature [^oC]')
        axis square; axis equal
        title({['Temperature evolution after ',...
            num2str(T.time(n)/T.year/1e6),' Myrs'];...
            ['for the ',FDSchema,' method and dt = ',...
            num2str(T.dt/(1/(2*P.kappa*(1/N.dx^2+1/N.dz^2)))),...
            '*dt_{critical}']})
        set(gca,'FontWeight','Bold','FontSize',10,'LineWidth',2)
        
        #         # Capture the plot as an image
        #         frame       = getframe(h);
        #         im          = frame2im(frame);
        #         [imind,cm]  = rgb2ind(im,256);
        #
        #         # Wri te to the GIF File
        #         if n == 1
        #             imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
        #         else
        #             imwrite(imind,cm,filename,'gif','WriteMode','append');
        #         end
    end
end
# tendt       = toc(tstartt); 
# ttotal      = cputime-tstartt; 
##

figure
clf
switch Model
    case 'Gaussian'        
        subplot(1,3,1)
        plot(D.TProfile(:,1:N.dn:end)-P.K0,N.z(:,1)./1e3,'k')
    otherwise
        subplot(1,2,1)
        plot(D.TProfile(:,1:N.dn:end)-P.K0,N.z(:,1)./1e3)
end

switch Model
    case 'Gaussian'
        hold on
        plot(D.TProfilea(:,1:N.dn:end)-P.K0,N.z(:,1)./1e3,'y--')
end
xlabel('T_{x = L/2} [ ^oC ]'); ylabel('Depth [ km ]')
title([{'Temperatur Profile'};{'@ Distance x = L/2'}])
set(gca,'FontWeight','Bold','FontSize',10,'LineWidth',2)

switch Model
    case 'Gaussian'        
        subplot(1,3,2)        
        plot(T.time/T.year/1e6,D.Tmax-P.K0,T.time/T.year/1e6,D.Tmaxa-P.K0)
    otherwise
        subplot(1,2,2)
        plot(T.time/T.year/1e6,D.Tmax-P.K0)
end
xlabel('Time [ My ]'); ylabel('T_{max} [ ^oC ]')
title([{'Maximum Temperature'};{'@ Depth z = H/2'}])
set(gca,'FontWeight','Bold','FontSize',10,'LineWidth',2)
switch Model
    case 'Gaussian'
        legend('num','ana')
end

switch Model
    case 'Gaussian'                
        subplot(2,3,3)
        pcolor(N.x/1e3,N.z/1e3,D.epsT); shading interp; colorbar        
        xlabel('x [km]'); ylabel('z [km]'); zlabel('Temperature [^oC]')
        axis square; axis equal
        title({['Temperature deviation after ',...
            num2str(T.time(n)/T.year/1e6),' Myrs'];...
            ['for the ',FDSchema,' method and dt = ',...
            num2str(T.dt/(1/(2*P.kappa*(1/N.dx^2+1/N.dz^2)))),...
            '*dt_{critical}']})
        set(gca,'FontWeight','Bold','FontSize',10,'LineWidth',2)
        axis([0 P.L/1e3 -P.H/1e3 0])
        
        subplot(2,3,6)
        plot(T.time/T.year/1e6,D.RMS)
        ylabel('RMS'); xlabel('Time[ My ]')        
        set(gca,'FontWeight','Bold','FontSize',10,'LineWidth',2)
end


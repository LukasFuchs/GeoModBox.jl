# Diffusionsgleichung (2D)
function Diffusion_2D_discretization_comparison
clear 
clc

Model       =   'Gaussian';
# 'explicit','implicit','CNV','ADI'
Schema    =   {'explicit';'implicit';'CNV';'ADI';'implicito'};

# Physikalischer Parameter ---------------------------------------------- #
P.L         =   200e3;          #   L??nge des Models    [m]
P.H         =   200e3;          #   H??he des Models     [m]

P.k         =   3;              #   Thermische Konduktivit??t [W/m/K]
P.cp        =   1000;           #   W??rmekapazit??t [J/kg/K]
P.rho       =   3200;           #   Dichte [kg/m^3]
P.kappa     =   P.k/P.rho/P.cp; #   Thermische Diffusivit??t [m^2/s]
P.K0        =   273.15;         #   Kelvin bei 0 C

P.Q0        =   0;              #   Hintergrund Waermeproduktionsrate;

# N.dn        =   5;              # Inkremente der Graphischen Darstellung, d.h. hier nur jeder 25 Zeitschritt

P.Tbot      =   500;            # Temperatur am unteren Rand  [C]
P.Tbot      =   P.Tbot+P.K0;    # Temperatur in Kelvin
P.Ttop      =   500;            # Temperatur am oberen Rand   [C]
P.Ttop      =   P.Ttop+P.K0;    # Temperatur in Kelvin

P.Tmax      =   500;            # Temperaturamplitude [C]
P.Tmax      =   P.Tmax+P.K0;    # Temperatur in Kelvin
P.sigma     =   8e3;            #
P.Xc        =   P.L/2;          # X-Koordinate des Zentrums
P.Zc        =   P.H/2;          # z-Koordinate des Zentrums

eps         = zeros(size(Schema,1),4);
epsr        = zeros(size(Schema,1),4);
nxnz        = zeros(size(Schema,1),4);
Tmax        = zeros(size(Schema,1),4);
Tmean       = zeros(size(Schema,1),4);

for k = 1:size(Schema,1)
    FDSchema = Schema{k};
    disp('')
    disp(FDSchema)
    for l = 1:3
        # Numerische Parameter ------------------------------------------ #
        N.nx        =   l*20+1;           # # Gitterpunkte in x-Richtung
        N.nz        =   l*20+1;           # # Gitterpunkte in z-Richtung
        N.dx        =   P.L/(N.nx-1);   # Gitterabstand in x-Richtung
        N.dz        =   P.H/(N.nz-1);   # Gitterabstand in z-Richtung
        
        disp(['nx = ',num2str(N.nx)])
        
        [N.x,N.z]   =   meshgrid(0:N.dx:P.L,-P.H:N.dz:0);
        
        N.beenhere  =   0;              # Dummy parameter, nicht ver??ndern!
        
        # Berechnung der stabilen Zeitschrittlaenge --------------------- #
        T.tmax      =   20;                 # Maximale Laufzeit des Models in Ma
        T.day       =   3600*24;            # # Sekunden pro Tag
        T.year      =   365.25*T.day;       # # Sekunden pro Jahr
        T.tmax      =   T.tmax*1e6*T.year;  # Maximale Zeit in Sekunden
        
        T.dtfac     =   1.0;                # Multiplikationsfaktor f??r dt
        T.dt        =   T.dtfac*(1/(2*P.kappa*(1/N.dx^2+1/N.dz^2)));
        
        N.nt        =   ceil(T.tmax/T.dt);  # Anzahl der Zeitschritte
        T.time      =   zeros(1,N.nt);
        #
        
        B.ttbc      =   'const';            # Oben - bisher nur const
        B.btbc      =   'const';            # Unten - bisher nur const
        B.ltbc      =   'const';             # links - const od. flux
        B.rtbc      =   'const';             # rechts - const od. flux
        ##
        
        # Waermefluss
        # ----------
        B.lhf       = 0;    B.rhf       = 0;
        # ----------
        B.bhf       = 0;    B.thf       = 0;
        ## Erstellung des Anfangstemperaturfseld
        
        D.T0    =  P.Tbot + P.Tmax .*...
            exp( -( (N.x-P.Xc).^2 + (N.z+P.Zc).^2 ) ./...
            ( 2*P.sigma^2/pi ));
        D.Tana          = D.T0;
        
        # Hintergrundfeld f??r W??rmequellen
        D.Q             = P.Q0.*ones(N.nz,N.nx);
                
        D.RMS           = zeros(1,N.nt);        
        
        ## Zeitschleife        
        for n=1:N.nt
            if n>1
                switch FDSchema
                    case 'explicit'
                        D.T1        =   ...
                            SolveDiff2Dexplicit(D.T0,D.Q,T.dt,P,N,B);
                        
                    case 'implicit'
                        [D.T1,N]    =   ...
                            SolveDiff2Dimplicit(D.T0,D.Q,T.dt,P,N,B);
                    case 'implicito'
                        [D.T1,N]    =   ...
                            SolveDiff2Dimplicit_opt(D.T0,D.Q,T.dt,P,N,B);                
                    case 'ADI'
                        if(n==1)
                            N.beenhere = 0;
                        end
                        [D.T1,N]    =   ...
                            SolveDiff2DADI(D.T0,D.Q,T.dt,P,N,B);
                    case 'CNV'
                        if(n==1)
                            N.beenhere = 0;
                        end
                        [D.T1,N]    = ...
                            SolveDiff2DCNV(D.T0,D.Q,T.dt,P,N,B);
                end
                
                # Zuweisung der neuen Temperatur
                D.T0        =   D.T1;
                
                T.time(n)   =   T.time(n-1) + T.dt;                
                
                D.Tana      =  P.Tbot + P.Tmax ./...
                    (1 + 2*pi*T.time(n)*P.kappa/P.sigma^2) .*...
                    exp( -( (N.x-P.Xc).^2 + (N.z+P.Zc).^2 ) ./...
                    ( 2*P.sigma^2/pi + 4*T.time(n)*P.kappa));
                
            end
            
            D.Tmax(n)           = max(max(D.T0));            
            D.Tmean(n)          = mean(mean(D.T0));            
            
            D.epsT      = (D.Tana-D.T0);
            D.epsTr     = (D.Tana-D.T0)./D.Tana;
            
            D.RMS(n)    = sqrt(sum(sum((D.epsT).^2))/(N.nx*N.nz));
            D.RMSr(n)   = sqrt(sum(sum((D.epsTr).^2))/(N.nx*N.nz));
            
            
            if (n==N.nt)
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
                
            end
        end
        
        disp(' Time loop finished ... ')
        disp('-> Use new grid size...')
        
        eps(k,l)    = D.RMS(N.nt);
        epsr(k,l)   = D.RMSr(N.nt);
        nxnz(k,l)   = N.nx*N.nz;
        Tmax(k,l)   = D.Tmax(N.nt);
        Tmean(k,l)  = D.Tmean(N.nt);
        
        Tana        = D.Tana; 
        clear N T B
    end
end

legendinfo  = cell(1,size(Schema,1)+1);
linstyle    = {'-','--',':','-.','-'};
figure(2)
clf
for k = 1:size(Schema,1)
    subplot(1,3,1)
    p(k) = loglog(1./nxnz(k,:),eps(k,:),'LineStyle',linstyle{k},...
        'Marker','*','LineWidth',2);
#     loglog(1./nxnz(k,:),epsr(k,:),'LineStyle',linstyle{k},...
#         'Marker','o','LineWidth',2);
    legendinfo{k} = Schema{k};
    hold on
    if k == size(Schema,1)        
        legendinfo{k+1} = '';
        xlabel('1/nx/nz'); ylabel('RMS_{\DeltaT}')
        axis square        
        set(gca,'FontWeight','Bold','FontSize',10,'LineWidth',2)
        legend(p,legendinfo,'Location','NorthWest')
        axis([1e-5 1e-2 1e-2 2])
    end
    subplot(1,3,2)
    loglog(1./nxnz(k,:),Tmax(k,:),'LineStyle',linstyle{k},...
        'Marker','*','LineWidth',2)
    hold on
    if k == size(Schema,1)
        legendinfo{k+1} = 'Sol_{ana}';
        plot(1/max(nxnz(1,:)):1e-4:1/min(nxnz(1,:)),...
            max(max(D.Tana)).*...
            ones(1,length(1/max(nxnz(1,:)):1e-4:1/min(nxnz(1,:)))),'k--')
        xlabel('1/nx/nz'); ylabel('T_{max}')
        axis square
        set(gca,'FontWeight','Bold','FontSize',10,'LineWidth',2,...
            'xscale','log')
        legend(legendinfo,'Location','NorthWest')
        axis([1e-5 1e-2 785 795])
    end
    subplot(1,3,3)
    loglog(1./nxnz(k,:),Tmean(k,:),'LineStyle',linstyle{k},...
        'Marker','*','LineWidth',2)
    hold on
    if k == size(Schema,1)
        plot(1/max(nxnz(1,:)):1e-4:1/min(nxnz(1,:)),...
            mean(mean(D.Tana)).*...
            ones(1,length(1/max(nxnz(1,:)):1e-4:1/min(nxnz(1,:)))),'k--')
        xlabel('1/nx/nz'); ylabel('\langle T \rangle')
        axis square
        set(gca,'FontWeight','Bold','FontSize',10,'LineWidth',2)
        legend(legendinfo,'Location','NorthEast')
        axis([1e-5 1e-2 775 776])
    end
end

# end


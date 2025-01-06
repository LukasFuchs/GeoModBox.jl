% function Avection2D
% ----------------------------------------------------------------------- %
% Funktion zur Loesung des zweidimensionalen Advektionsproblem mit Hilf von
% unterschiedlichen Methoden. Zu waehlen sind:
%   'upwind'    - Upwind Schema
%   'slf'       - Staggered Leaped Frog Schema
%   'semi-lag'  - Semi-Lagrangian Schema
%   'tracers'   - Tracer Methode (vollkommen Lagrangian)
%
% Für die Temperaturanomalie kann entweder ein rechteckiger Block gewählt
% werden oder eine Gausssche Temperaturverteilung:
%   'block'     - Rechteckiger Block
%   'gaussian'  - Gaussche Temperaturverteilung
%   'circle'    - Kreisförmige Anomalie
%
% Für das konstante Geschwindigkeitsfeld könne  n zwei Varianten gewählt
% werden:
%   'RigidBody'     - Ein Rotationsfeld mit konstanter Rotation
%   'ShearCell'     - Eine konstantes Konvektionsfeld mit Sheardeformation
%
% ----------------------------------------------------------------------- %
% Vers. 1.0 - 1.04.2020
% ======================================================================= %
clear
clc
profile on

fdmethod    = 'upwind';
Tanomaly    = 'block';
FlowField   = 'RigidBody';

%% Plot Settings -------------------------------------------------------- %
inc     = 10; % Increment for quiver vectors

%% Geometrical constants ------------------------------------------------ %
xmin    = 0;
xmax    = 1;
zmax    = xmax;
zmin    = xmin;

%% Numerical constants -------------------------------------------------- %
nx      = 101;
nz      = nx;
dx      = (abs(xmin)+xmax)/(nx-1);
dz      = (abs(zmin)+zmax)/(nz-1);

x       = xmin:dx:xmax;
z       = xmin:dz:zmax;

[X,Z]   = meshgrid(x,z);

%% Zeit Parameter ------------------------------------------------------- %
tmax    = 6.336;

dtfac   = 1.0;           % Courant time factor, i.e. dtfac*dt_courant

%% Tracer advection method ---------------------------------------------- %
nmx     = 5;
nmz     = 5;

%% Festlegung der Temperaturanomalie ------------------------------------ %
switch Tanomaly
    case 'circle'
        % Hintergrund Temperatur ---------------------------------------- %
        Tb      =   1000;
        
        % Anomalie Temperatur ------------------------------------------- %
        Ta      =   1200;
        T       =   zeros(nz,nx);
        
        % Bereich der Anomalie ------------------------------------------ %
        ri          =   .2;
        xc          =   (xmin+xmax)/4;
        zc          =   (zmin+zmax)/2;
        alpha       =   0.0;
        a_ell       =   .2;
        b_ell       =   .2;
        x_ell       =   X .* cosd(alpha) + Z .* sind(alpha);
        z_ell       =   -X .* sind(alpha) + Z .* cosd(alpha);
        Elli        =   ((x_ell - xc)./ a_ell).^2 + ((z_ell-zc)./ b_ell).^2;
        T(Elli<=ri) =   Ta;
        T(Elli>ri)  =   Tb;
        
        Tmax        =   max(max(T));
        Tmin        =   min(min(T));
        Tmean       =   (Tmax+Tmin)/2;
    case 'gaussian'
        % Gaussche Temperatur Anomalie ---------------------------------- %
        Ampl    = 200;     % Amplitude der Anomalie
        sigma   = 0.1;      % Breite der Anomalie
        T0      = 1000;     % Hintergrund Temperatur
        
        T       = T0 + Ampl*exp(-((X - 0.25).^2 + (Z - 0.5).^2)./sigma^2);
        
        Tmax    = max(max(T));
        Tmin    = min(min(T));
        Tmean   = (Tmax+Tmin)/2;
        
    case 'block'
        % Hintergrund Temperatur ---------------------------------------- %
        Tb      = 1000;
        
        % Bereich der Temperatur Anomalie ------------------------------- %
        xTl     = (abs(xmin)+abs(xmax))/4 - (abs(xmin)+abs(xmax))/10;
        xTr     = (abs(xmin)+abs(xmax))/4 + (abs(xmin)+abs(xmax))/10;
        zTu     = (abs(zmin)+abs(zmax))/2 - (abs(zmin)+abs(zmax))/10;
        zTo     = (abs(zmin)+abs(zmax))/2 + (abs(zmin)+abs(zmax))/10;
        
        Ta      = 1200;
        
        Tmean   = (Tb + Ta)/2;
        
        % Anfangstemperatur Verteilung ---------------------------------- %
        T       = Tb.*ones(nz,nx);
        T(z>=zTu&z<=zTo,x>=xTl&x<=xTr) = Ta;
        %         tc      = 100;
        Tmax    = max(max(T));
end

switch fdmethod
    case 'slf'
        Told    = T;
    case 'tracers'
        nmxx        = (nx-1)*nmx;
        nmzz        = (nz-1)*nmz;
        dmx         = (abs(xmin)+abs(xmax))/(nmxx-1);
        dmz         = (abs(zmin)+abs(zmax))/(nmzz-1);
        xm          = linspace(xmin,xmax-dmx,nmxx);
        zm          = linspace(zmin,zmax-dmz,nmzz);
        [XM,ZM]     = meshgrid(xm,zm);
        XM          = XM + rand(nmzz,nmxx)*dmx;
        XM          = reshape(XM,[nmzz*nmxx,1]);
        ZM          = ZM + rand(nmzz,nmxx)*dmz;
        ZM          = reshape(ZM,[nmzz*nmxx,1]);
        Tm          = zeros(nmzz,nmxx);
        Tm          = reshape(Tm,[nmzz*nmxx,1]);
        [Tm,~]      = TracerInterp(Tm,XM,ZM,T,[],X,Z,'to');
end

%% Definition des Geschwindigkeitsfeldes -------------------------------- %
switch FlowField
    case 'RigidBody'
        % Starre Rotation
        vx      = (Z-(zmax-zmin)/2);            % [ m/s ]
        vz      = -(X-(xmax-xmin)/2);           % [ m/s ]
        
        switch fdmethod
            case 'tracers'
                Rad = sqrt((X-(xmax-xmin)/2).^2 + (Z-(zmax-zmin)/2).^2);
                
                vx(Rad>((xmax-xmin)/2)) = 0;
                vz(Rad>((xmax-xmin)/2)) = 0;
        end
    case 'ShearCell'
        % Zelle mit einfacher Scherdehnung
        vx      = -sin(pi.*X).*cos(pi.*Z);
        vz      = cos(pi.*X).*sin(pi.*Z);
end

% Define time step ------------------------------------------------------ %
% dt      = dtfac*min(dx,dz)/max(max(max(abs(vx))),max(max(abs(vz))));
dt      = dtfac*min(dx,dz)/(sqrt(max(max(vx))^2 + max(max(vz))^2));
nt      = ceil(tmax/dt);

%% Animation settings
filename    = ['2D_Advection_',fdmethod,'_',Tanomaly,'_',FlowField,'.gif'];
if strcmp(fdmethod,'tracers')
    if strcmp(getenv('OS'),'Windows_NT')
        set(figure(1),'position',[55.4,125.8,1326.4,636.2]);
    else
        set(figure(1),'position',[170,221,1558,684]);
    end
end
h           = figure(1);

% Solve advection equation ---------------------------------------------- %
for t=1:nt
    disp([' Time step: ',num2str(t)])
    
    if t > 1
        switch fdmethod
            case 'upwind'
                T       = UpwindAdvection2D(vx,vz,T,dx,dz,dt);
            case 'slf'
                Tnew    = SLFAdvection2D(vx,vz,Told,T,dx,dz,dt);
                Told    = T;
                
                T       = Tnew;
            case 'semi-lag'
                Tnew    = SemiLagAdvection2D(vx,vz,[],[],X,Z,T,dt);
                
                T       = Tnew;
            case 'tracers'
                % if (t>2)
                %     % Interpolate temperature grid to tracers --------------- %
                %     [Tm,~]  = TracerInterp(Tm,XM,ZM,T,[],X,Z,'to');
                % end
                
                % Advect tracers with Runge-Kutta 4th order ----------------- %
                [XM,ZM] = ...
                    AdvectMarker2D(X,Z,XM,ZM,dt,vx,vz,xmax,xmin,zmax,zmin);
                
                % Interpolate temperature from tracers to grid -------------- %
                Told    = T;
                [~,T]   = TracerInterp(Tm,XM,ZM,T,[],X,Z,'from');
        end
    end
    
    if (mod(t,5)==0||t==1)
        switch fdmethod
            case 'tracers'
                figure(1),clf
                subplot(1,2,1)
                pcolor(X,Z,T./Tmax)
                shading interp; lighting phong; hold on;
                k = colorbar('southoutside');
                title(k,'T / T_{max}','Position',[320 0])
                quiver(X(1:inc:end,1:inc:end),Z(1:inc:end,1:inc:end),...
                    vx(1:inc:end,1:inc:end),vz(1:inc:end,1:inc:end))
                xlabel('x [ m ]'); ylabel('z [ m ]')
                axis equal; axis tight
                title({['2-D numerical Advection:',fdmethod];...
                    ['\Deltat_{fac} = ',num2str(dtfac),...
                    '; nx = ',num2str(nx),', nz = ',num2str(nz),', mpe: ',num2str(nmx)];...
                    ['Step: ',num2str(t)]})
                
                subplot(1,2,2)
                plot(XM,ZM,'.','MarkerSize',1)
                hold on
                plot(X,Z,'kx','MarkerSize',2)
                contour(X,Z,T,[Tmean Tmean],'k','LineWidth',1)
                xlabel('x [ m ]'); ylabel('z [ m ]')
                title('Tracerdistribution')
                axis equal; axis tight
                drawnow
            otherwise
                figure(1),clf
                pcolor(X,Z,T./Tmax)
                shading interp; lighting phong; hold on
                k = colorbar;
                title(k,'T / T_{max}')
                quiver(X(3:inc:end-2,3:inc:end-2),Z(3:inc:end-2,3:inc:end-2),...
                    vx(3:inc:end-2,3:inc:end-2),vz(3:inc:end-2,3:inc:end-2))
                xlabel('x [ m ]'); ylabel('z [ m ]')
                title({['2-D numerical Advection:',fdmethod];...
                    ['\Deltat_{fac} = ',num2str(dtfac),...
                    '; nx = ',num2str(nx),', nz = ',num2str(nz)];...
                    ['Step: ',num2str(t)]})
                axis equal; axis tight
                drawnow
        end
        
        %         % Capture the plot as an image
        %         frame       = getframe(h);
        %         im          = frame2im(frame);
        %         [imind,cm]  = rgb2ind(im,256);
        %
        %         % Write to the GIF File
        %         if t == 1
        %             imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
        %         else
        %             imwrite(imind,cm,filename,'gif','WriteMode','append');
        %         end
    end
end
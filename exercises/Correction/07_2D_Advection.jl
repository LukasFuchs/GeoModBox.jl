# -------------------------------------------------------------------- #
# Funktion zur Loesung des zweidimensionalen Advektionsproblem mit 
# Hilfe von unterschiedlichen Methoden. Zu waehlen sind:
#   'upwind'    - Upwind Schema
#   'slf'       - Staggered Leaped Frog Schema
#   'semi-lag'  - Semi-Lagrangian Schema
#   'tracers'   - Tracer Methode (vollkommen Lagrangian)
#
# Für die Temperaturanomalie kann entweder ein rechteckiger Block 
# gewählt werden oder eine Gausssche Temperaturverteilung:
#   'block'     - Rechteckiger Block
#   'gaussian'  - Gaussche Temperaturverteilung
#   'circle'    - Kreisförmige Anomalie
#
# Für das konstante Geschwindigkeitsfeld könne  n zwei Varianten 
# gewählt werden:
#   'RigidBody' - Ein Rotationsfeld mit konstanter Rotation
#   'ShearCell' - Eine konstantes Konvektionsfeld mit Sheardeformation
#
# -------------------------------------------------------------------- #
# Vers. 1.0 - 26.11.2024 - Julia
# ==================================================================== #
# using Statistics
using Plots

function Advection_2D()

# Definition numerischer Verfahren =================================== #
FD          =   (Method     = (Adv=:upwind,),)
Ini         =   (T=:circle,V=:RigidBody,)
# -------------------------------------------------------------------- #

#fdmethod    =:upwind
#Tanomaly    =:block
#FlowField   =:RigidBody

# Plot Einstellungen ================================================= #
Pl  =   (
    inc         =   10,
)
# -------------------------------------------------------------------- #
# Model Konstanten =================================================== #
M   =   (
    xmin    =   0.0,
    xmax    =   1.0,
    ymin    =   0.0,
    ymax    =   1.0,
)
# -------------------------------------------------------------------- #
# Numerische Konstanten ============================================== #
NC  =   (
    xc      =   100,        # Number of horizontal centroids
    yc      =   100,        # Number of vertical centroids
)
NC1 =   (
    xv      =   NC.xc + 1,  # Number of horizontal vertices
    yv      =   NC.yc + 1,  # Number of vertical vertices
)
NC  =   merge(NC,NC1)

Δ   =   (
    x   =   (abs(M.xmin)+M.xmax)/NC.xc,
    y   =   (abs(M.ymin)+M.ymax)/NC.yc,
)
# -------------------------------------------------------------------- #
# Erstellung des Gitters ============================================= #
x   =   (
    c       =   LinRange(M.xmin + Δ.x/2.0, M.xmax - Δ.x/2.0, NC.xc),
    v       =   LinRange(M.xmin, M.xmax , NC.xv)
)
y       = (
    c       =   LinRange(M.ymin + Δ.y/2.0, M.ymax - Δ.y/2.0, NC.yc),
    v       =   LinRange(M.ymin, M.ymax, NC.yv),
)
xc2d = x.c .+ 0*y.c'
yc2d = 0*x.c .+ y.c'
# -------------------------------------------------------------------- #
# Zeit Konstanten ==================================================== #
T   =   ( 
    tmax    = 6.336,
    dtfac   = 1.0,          # Courant time factor, i.e. dtfac*dt_courant
)
# -------------------------------------------------------------------- #
# Tracer Advektions Verfahren ======================================== #
Ma  =   (
    nmx     =   5,
    nmz     =   5
)
# -------------------------------------------------------------------- #
# Anfangsbedingungen ================================================= #
# Temperatur --------------------------------------------------------- #
D       =   (
    T       =   zeros(NC.xc,NC.yc),
    vx      =   zeros(NC.xv+1,NC.yv),
    vy      =   zeros(NC.xv,NC.yv+1),
    vxc     =   zeros(NC.xc,NC.yc),
    vyc     =   zeros(NC.xc,NC.yc),
    vc      =   zeros(NC.xc,NC.yc),
    Tmax    =   [0.0],
    Tmin    =   [0.0],
    Tmean   =   [0.0],
)
if Ini.T==:circle    
    # Hintergrund Temperatur ----------------------------------------- #
    Tb          =   1000
        
    # Anomalie Temperatur -------------------------------------------- #
    Ta          =   1200
        
    # Bereich der Anomalie ------------------------------------------- #
    ri          =   .2
    xc          =   (M.xmin+M.xmax)/4
    yc          =   (M.ymin+M.ymax)/2
    α           =   0.0
    a_ell       =   .2
    b_ell       =   .2
    for i = 1:NC.xc, j = 1:NC.yc
        x_ell   =  x.c[i]*cosd(α) + y.c[j]*sind(α)
        y_ell   =  -x.c[i]*sind(α) + y.c[j]*cosd(α)
        Elli    =   ((x_ell - xc)/ a_ell)^2 + ((y_ell-yc)/ b_ell)^2
        if Elli <= ri 
            D.T[i,j]    =   Ta
        else
            D.T[i,j]    =   Tb
        end
    end
        
    D.Tmax[1]   =   maximum(D.T)
    D.Tmin[1]   =   minimum(D.T)
    D.Tmean[1]  =   (D.Tmax[1]+D.Tmin[1])/2
elseif Ini.T==:gaussian
    # Gaussche Temperatur Anomalie ----------------------------------- #
    Ampl        =   200     # Amplitude der Anomalie
    σ           =   0.1     # Breite der Anomalie
    T0          =   1000    # Hintergrund Temperatur
    κ           =   1e-6
        
    # AnalyticalSolution2D!(D.T, x.c, y.c, 0.0, (T0=Ampl,K=κ,σ=σ))

    for i = 1:NC.xc, j = 1:NC.yc
        D.T[i,j]    =   T0 + Ampl*exp(-((x.c[i] - 0.25)^2 + (y.c[j] - 0.5)^2)/σ^2)
    end
        
    D.Tmax[1]   =   maximum(D.T)
    D.Tmin[1]   =   minimum(D.T)
    D.Tmean[1]  =   (D.Tmax[1]+D.Tmin[1])/2
        
elseif Ini.T==:block
    # Hintergrund Temperatur ----------------------------------------- #
    Tb      =   1000
        
    # Bereich der Temperatur Anomalie -------------------------------- #
    xTl     =   (abs(M.xmin)+abs(M.xmax))/4 - (abs(M.xmin)+abs(M.xmax))/10
    xTr     =   (abs(M.xmin)+abs(M.xmax))/4 + (abs(M.xmin)+abs(M.xmax))/10
    yTu     =   (abs(M.ymin)+abs(M.ymax))/2 - (abs(M.ymin)+abs(M.ymax))/10
    yTo     =   (abs(M.ymin)+abs(M.ymax))/2 + (abs(M.ymin)+abs(M.ymax))/10
        
    Ta      =   1200
        
    D.Tmean[1]  =   (Tb + Ta)/2
         
    # Anfangstemperatur Verteilung ----------------------------------- #
    for i = i:NC.xc, j = 1:NC.yc
        if y.c[j]>=yTu & y.c[j] <= yTo & x.c[i]>=xTl & x.c[i]<=xTr
            D.T[i,j]    =   Ta
        end
    end
    D.Tmax[1]    = maximum(D.T)
end

# if FD.Method.Adv==slf
#     Told    = T;
# elseif FD.Method.Adv==tracers
#     nmxx        = (nx-1)*nmx;
#     nmzz        = (nz-1)*nmz;
#     dmx         = (abs(xmin)+abs(xmax))/(nmxx-1);
#     dmz         = (abs(zmin)+abs(zmax))/(nmzz-1);
#     xm          = linspace(xmin,xmax-dmx,nmxx);
#     zm          = linspace(zmin,zmax-dmz,nmzz);
#     [XM,ZM]     = meshgrid(xm,zm);
#     XM          = XM + rand(nmzz,nmxx)*dmx;
#     XM          = reshape(XM,[nmzz*nmxx,1]);
#     ZM          = ZM + rand(nmzz,nmxx)*dmz;
#     ZM          = reshape(ZM,[nmzz*nmxx,1]);
#     Tm          = zeros(nmzz,nmxx);
#     Tm          = reshape(Tm,[nmzz*nmxx,1]);
#     [Tm,~]      = TracerInterp(Tm,XM,ZM,T,[],X,Z,'to');
# end

# Definition des Geschwindigkeitsfeldes ------------------------------ #
if Ini.V==:RigidBody
    # Starre Rotation
    for i = 1:NC.xv+1, j = 1:NC.yv
        D.vx[i,j]  =    (y.v[j]-(M.ymax-M.ymin)/2)            # [ m/s ]
    end
    for i = 1:NC.xv, j = 1:NC.yv+1
        D.vy[i,j]  =   -(x.v[i]-(M.xmax-M.xmin)/2)           # [ m/s ]
    end
        
    #if FD.Method.Adv==:tracers
    #        Rad = sqrt((X-(xmax-xmin)/2).^2 + (Z-(zmax-zmin)/2).^2);
    #            
    #        vx(Rad>((xmax-xmin)/2)) = 0;
    #        vz(Rad>((xmax-xmin)/2)) = 0;
    #end
elseif Ini.V==:ShearCell
    # Zelle mit einfacher Scherdehnung
    for i = 2:NC.xv j = 1:NC.yv
        D.vx[i,j]   =   -sin(π*x.c[i-1])*cos(π*y.v[j])
    end
    for i = 1:NC.xv, j = 1:NC.yv+1
        D.vz[i,j]   =   cos(pi.*X).*sin(pi.*Z);
    end
end
for i = 1:NC.xc, j = 1:NC.yc
    D.vxc[i,j]  = (D.vx[i+1,j] + D.vx[i+1,j+1])/2
    D.vyc[i,j]  = (D.vy[i,j+1] + D.vy[i+1,j+1])/2
end
@. D.vc        = sqrt(D.vxc^2 + D.vyc^2)

# Visualize initial condition ---
p = heatmap(x.c , y.c, D.T', 
        color=:thermal, colorbar=false, aspect_ratio=:equal, 
        xlabel="x", ylabel="z", 
        title="Temperature", 
        xlims=(M.xmin, M.xmax), ylims=(M.ymin, M.ymax), 
        clims=(minimum(D.T), maximum(D.T)),layout=(1,2),
        subplot=1)
# quiver!(p,xc2d,yc2d,quiver=(D.vxc,D.vyc),subplot=1)
heatmap!(p,x.c,y.c,D.vc',
        color=:viridis, colorbar=true, aspect_ratio=:equal, 
        xlims=(M.xmin, M.xmax), ylims=(M.ymin, M.ymax),
        subplot=2)
# -------------------------------------------------------------------- #
## Define time step ------------------------------------------------------ #
## dt      = dtfac*min(dx,dz)/max(max(max(abs(vx))),max(max(abs(vz))));
#dt      = dtfac*min(dx,dz)/(sqrt(max(max(vx))^2 + max(max(vz))^2));
#nt      = ceil(tmax/dt);

# ## Animation settings
# filename    = ['2D_Advection_',fdmethod,'_',Tanomaly,'_',FlowField,'.gif'];
# if strcmp(fdmethod,'tracers')
#     if strcmp(getenv('OS'),'Windows_NT')
#         set(figure(1),'position',[55.4,125.8,1326.4,636.2]);
#     else
#         set(figure(1),'position',[170,221,1558,684]);
#     end
# end
# h           = figure(1);

# # Solve advection equation ---------------------------------------------- #
# for t=1:nt
#     disp([' Time step: ',num2str(t)])
    
#     if t > 1
#         switch fdmethod
#             case 'upwind'
#                 T       = UpwindAdvection2D(vx,vz,T,dx,dz,dt);
#             case 'slf'
#                 Tnew    = SLFAdvection2D(vx,vz,Told,T,dx,dz,dt);
#                 Told    = T;
                
#                 T       = Tnew;
#             case 'semi-lag'
#                 Tnew    = SemiLagAdvection2D(vx,vz,[],[],X,Z,T,dt);
                
#                 T       = Tnew;
#             case 'tracers'
#                 # if (t>2)
#                 #     # Interpolate temperature grid to tracers --------------- #
#                 #     [Tm,~]  = TracerInterp(Tm,XM,ZM,T,[],X,Z,'to');
#                 # end
                
#                 # Advect tracers with Runge-Kutta 4th order ----------------- #
#                 [XM,ZM] = ...
#                     AdvectMarker2D(X,Z,XM,ZM,dt,vx,vz,xmax,xmin,zmax,zmin);
                
#                 # Interpolate temperature from tracers to grid -------------- #
#                 Told    = T;
#                 [~,T]   = TracerInterp(Tm,XM,ZM,T,[],X,Z,'from');
#         end
#     end
    
#     if (mod(t,5)==0||t==1)
#         switch fdmethod
#             case 'tracers'
#                 figure(1),clf
#                 subplot(1,2,1)
#                 pcolor(X,Z,T./Tmax)
#                 shading interp; lighting phong; hold on;
#                 k = colorbar('southoutside');
#                 title(k,'T / T_{max}','Position',[320 0])
#                 quiver(X(1:inc:end,1:inc:end),Z(1:inc:end,1:inc:end),...
#                     vx(1:inc:end,1:inc:end),vz(1:inc:end,1:inc:end))
#                 xlabel('x [ m ]'); ylabel('z [ m ]')
#                 axis equal; axis tight
#                 title({['2-D numerical Advection:',fdmethod];...
#                     ['\Deltat_{fac} = ',num2str(dtfac),...
#                     '; nx = ',num2str(nx),', nz = ',num2str(nz),', mpe: ',num2str(nmx)];...
#                     ['Step: ',num2str(t)]})
                
#                 subplot(1,2,2)
#                 plot(XM,ZM,'.','MarkerSize',1)
#                 hold on
#                 plot(X,Z,'kx','MarkerSize',2)
#                 contour(X,Z,T,[Tmean Tmean],'k','LineWidth',1)
#                 xlabel('x [ m ]'); ylabel('z [ m ]')
#                 title('Tracerdistribution')
#                 axis equal; axis tight
#                 drawnow
#             otherwise
#                 figure(1),clf
#                 pcolor(X,Z,T./Tmax)
#                 shading interp; lighting phong; hold on
#                 k = colorbar;
#                 title(k,'T / T_{max}')
#                 quiver(X(3:inc:end-2,3:inc:end-2),Z(3:inc:end-2,3:inc:end-2),...
#                     vx(3:inc:end-2,3:inc:end-2),vz(3:inc:end-2,3:inc:end-2))
#                 xlabel('x [ m ]'); ylabel('z [ m ]')
#                 title({['2-D numerical Advection:',fdmethod];...
#                     ['\Deltat_{fac} = ',num2str(dtfac),...
#                     '; nx = ',num2str(nx),', nz = ',num2str(nz)];...
#                     ['Step: ',num2str(t)]})
#                 axis equal; axis tight
#                 drawnow
#         end
        
#         #         # Capture the plot as an image
#         #         frame       = getframe(h);
#         #         im          = frame2im(frame);
#         #         [imind,cm]  = rgb2ind(im,256);
#         #
#         #         # Write to the GIF File
#         #         if t == 1
#         #             imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
#         #         else
#         #             imwrite(imind,cm,filename,'gif','WriteMode','append');
#         #         end
#     end
# end

end # Function end

Advection_2D()
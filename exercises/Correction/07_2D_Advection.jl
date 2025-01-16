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
using GeoModBox.AdvectionEquation.TwoD
using GeoModBox.InitialCondition

function Advection_2D()

# Definition numerischer Verfahren =================================== #
# Define Advection Scheme ---
#   1) upwind, 2) slf
FD          =   (Method     = (Adv=:slf,),)
# Define Initial Condition ---
# Temperature - 
#   1) circle, 2) gaussian, 3) block
# Velocity - 
#   1) RigidBody, 2) ShearCell
Ini         =   (T=:circle,V=:RigidBody,) 
# -------------------------------------------------------------------- #
# Plot Einstellungen ================================================= #
Pl  =   (
    inc         =   5,
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
    cew     =   LinRange(M.xmin - Δ.x/2.0, M.xmax + Δ.x/2.0, NC.xc+2),
    v       =   LinRange(M.xmin, M.xmax , NC.xv)
)
y       = (
    c       =   LinRange(M.ymin + Δ.y/2.0, M.ymax - Δ.y/2.0, NC.yc),
    cns     =   LinRange(M.ymin - Δ.x/2.0, M.ymax + Δ.x/2.0, NC.yc+2),
    v       =   LinRange(M.ymin, M.ymax, NC.yv),
)
x1      =   ( 
    xc2d    =   x.c .+ 0*y.c',
    xv2d    =   x.v .+ 0*y.v', 
    xcew2d  =   x.cew .+ 0*y.v',
)
x   =   merge(x,x1)
y1      =   (
    yc2d    =   0*x.c .+ y.c',
    yv2d    =   0*x.v .+ y.v',
    ycns2d  =   0*x.v .+ y.cns',
)
y   =   merge(y,y1)
# -------------------------------------------------------------------- #
# Tracer Advektions Verfahren ======================================== #
Ma  =   (
    nmx     =   5,
    nmz     =   5
)
# -------------------------------------------------------------------- #
# Animationssettings ================================================= #
path        =   string("./exercises/Correction/Results/")
anim        =   Plots.Animation(path, String[] )
filename    =   string("07_2D_advection_",Ini.T,"_",Ini.V,
                        "_",FD.Method.Adv)
save_fig    =   1
# -------------------------------------------------------------------- #
# Anfangsbedingungen ================================================= #
# Temperatur --------------------------------------------------------- #
D       =   (
    T       =   zeros(NC.xc,NC.yc),
    T_ext   =   zeros(NC.xc+2,NC.yc+2),
    T_exto  =   zeros(NC.xc+2,NC.yc+2),
    vx      =   zeros(NC.xv,NC.yv+1),
    vy      =   zeros(NC.xv+1,NC.yv),
    vxc     =   zeros(NC.xc,NC.yc),
    vyc     =   zeros(NC.xc,NC.yc),
    vc      =   zeros(NC.xc,NC.yc),
    Tmax    =   [0.0],
    Tmin    =   [0.0],
    Tmean   =   [0.0],
)
if Ini.T==:circle    
    # Hintergrund Temperatur ---
    Tb          =   1000
    # Anomalie Temperatur ---
    Ta          =   1200
    # Bereich der Anomalie ---
    ri          =   .2
    xc          =   (M.xmin+M.xmax)/4
    yc          =   (M.ymin+M.ymax)/2
    α           =   0.0
    a_ell       =   .2
    b_ell       =   .2
    for i = 1:NC.xc+2, j = 1:NC.yc+2
        x_ell   =  x.cew[i]*cosd(α) + y.cns[j]*sind(α)
        y_ell   =  -x.cew[i]*sind(α) + y.cns[j]*cosd(α)
        Elli    =   ((x_ell - xc)/ a_ell)^2 + ((y_ell-yc)/ b_ell)^2
        if Elli <= ri 
            D.T_ext[i,j]    =   Ta
        else
            D.T_ext[i,j]    =   Tb
        end
    end    
    D.Tmax[1]   =   maximum(D.T_ext)
    D.Tmin[1]   =   minimum(D.T_ext)
    D.Tmean[1]  =   (D.Tmax[1]+D.Tmin[1])/2
elseif Ini.T==:gaussian
    # Gaussche Temperatur Anomalie ---
    Ampl        =   200     # Amplitude der Anomalie
    σ           =   0.05     # Breite der Anomalie
    T0          =   1000    # Hintergrund Temperatur
    # κ           =   1e-6
    # AnalyticalSolution2D!(D.T, x.c, y.c, 0.0, (T0=Ampl,K=κ,σ=σ))
    for i = 1:NC.xc+2, j = 1:NC.yc+2
        D.T_ext[i,j]    =   T0 + Ampl*exp(-((x.cew[i] - 0.20)^2 + (y.cns[j] - 0.5)^2)/σ^2)
    end
    # D.T         .=  D.T_ext[2:end-1,2:end-1]
    D.Tmax[1]   =   maximum(D.T_ext)
    D.Tmin[1]   =   minimum(D.T_ext)
    D.Tmean[1]  =   (D.Tmax[1]+D.Tmin[1])/2
elseif Ini.T==:block
    # Hintergrund Temperatur ---
    Tb      =   1000
    # Bereich der Temperatur Anomalie ---
    xTl     =   (abs(M.xmin-Δ.x/2)+abs(M.xmax+Δ.x/2))/4 - (abs(M.xmin-Δ.x/2)+abs(M.xmax+Δ.x/2))/10
    xTr     =   (abs(M.xmin-Δ.x/2)+abs(M.xmax+Δ.x/2))/4 + (abs(M.xmin-Δ.x/2)+abs(M.xmax+Δ.x/2))/10
    yTu     =   (abs(M.ymin-Δ.y/2)+abs(M.ymax+Δ.y/2))/2 - (abs(M.ymin-Δ.y/2)+abs(M.ymax+Δ.y/2))/10
    yTo     =   (abs(M.ymin-Δ.y/2)+abs(M.ymax+Δ.y/2))/2 + (abs(M.ymin-Δ.y/2)+abs(M.ymax+Δ.y/2))/10
    Ta      =   1200
    D.Tmean[1]  =   (Tb + Ta)/2
    # Anfangstemperatur Verteilung ---
    for i = 1:NC.xc+2, j = 1:NC.yc+2
        if y.cns[j]>=yTu && y.cns[j] <= yTo && x.cew[i]>=xTl && x.cew[i]<=xTr
            D.T_ext[i,j]    =   Ta
        end
    end
    #D.T         .=  D.T_ext[2:end-1,2:end-1]
    D.Tmax[1]   =   maximum(D.T_ext)
end
D.T         .=  D.T_ext[2:end-1,2:end-1]
if FD.Method.Adv==:slf
    D.T_exto    .=  D.T_ext
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
# Geschwindigkeit ---------------------------------------------------- #
IniVelocity!(Ini.V,D,NC,M,x,y)
# Get the velocity on the centroids ---
for i = 1:NC.xc, j = 1:NC.yc
    D.vxc[i,j]  = (D.vx[i,j+1] + D.vx[i+1,j+1])/2
    D.vyc[i,j]  = (D.vy[i+1,j] + D.vy[i+1,j+1])/2
end
@. D.vc        = sqrt(D.vxc^2 + D.vyc^2)
# Visualize initial condition ---------------------------------------- #
p = heatmap(x.c , y.c, (D.T./D.Tmax)', 
        color=:thermal, colorbar=true, aspect_ratio=:equal, 
        xlabel="x", ylabel="z", 
        title="Temperature", 
        xlims=(M.xmin, M.xmax), ylims=(M.ymin, M.ymax), 
        clims=(minimum(D.T)./D.Tmax[1], 1.0))
quiver!(p,x.xc2d[1:Pl.inc:end,1:Pl.inc:end],y.yc2d[1:Pl.inc:end,1:Pl.inc:end],
        quiver=(D.vxc[1:Pl.inc:end,1:Pl.inc:end],D.vyc[1:Pl.inc:end,1:Pl.inc:end]),        
        color="white")
if save_fig == 1
    Plots.frame(anim)
elseif save_fig == 0
    display(p)
end
# -------------------------------------------------------------------- #
# Zeit Konstanten ==================================================== #
T   =   ( 
    tmax    =   6.336,
    Δfac    =   1.0,    # Courant time factor, i.e. dtfac*dt_courant
    Δ       =   [0.0],
)
T.Δ[1]  =   T.Δfac * minimum((Δ.x,Δ.y)) / 
        (sqrt(maximum(abs.(D.vx))^2 + maximum(abs.(D.vy))^2))
nt      =   ceil(Int,T.tmax/T.Δ[1])
# -------------------------------------------------------------------- #
# Solve advection equation ------------------------------------------- #
 for i=2:nt
    display(string("Time step: ",i))    

    if FD.Method.Adv==:upwind
        upwindc2D!(D,NC,T,Δ)
    elseif FD.Method.Adv==:slf
        slfc2D!(D,NC,T,Δ)   
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
    end

    display(string("ΔT = ",((D.Tmax[1]-maximum(D.T))/D.Tmax[1])*100))

    # Plot Solution ---
    if mod(i,5) == 0 || i == nt
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

        p = heatmap(x.c , y.c, (D.T./D.Tmax)', 
                color=:thermal, colorbar=true, aspect_ratio=:equal, 
                xlabel="x", ylabel="z", 
                title="Temperature", 
                xlims=(M.xmin, M.xmax), ylims=(M.ymin, M.ymax), 
                clims=(minimum(D.T)./D.Tmax[1], 1.0))
        quiver!(p,x.xc2d[1:Pl.inc:end,1:Pl.inc:end],
                    y.yc2d[1:Pl.inc:end,1:Pl.inc:end],
                    quiver=(D.vxc[1:Pl.inc:end,1:Pl.inc:end].*0.2,
                            D.vyc[1:Pl.inc:end,1:Pl.inc:end].*0.2),        
                color="white")
        if save_fig == 1
            Plots.frame(anim)
        elseif save_fig == 0
            display(p)                        
        end
    end
        
end
# Save Animation ----------------------------------------------------- #
if save_fig == 1
    # Write the frames to a GIF file
    Plots.gif(anim, string( path, filename, ".gif" ), fps = 15)
    foreach(rm, filter(startswith(string(path,"00")), readdir(path,join=true)))
elseif save_fig == 0
    display(plot(p))
end

end # Function end

Advection_2D()
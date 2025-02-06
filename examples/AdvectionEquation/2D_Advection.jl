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
using Plots, Interpolations
using GeoModBox.AdvectionEquation.TwoD, GeoModBox.Tracers.TwoD
using GeoModBox.InitialCondition

function Advection_2D()

# Definition numerischer Verfahren =================================== #
# Define Advection Scheme ---
#   1) upwind, 2) slf, 3) semilag, 4) tracers
FD          =   (Method     = (Adv=:semilag,),)
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
    sc          =   1.0e9,
    Minc        =   1, 
    Msz         =   0.2,
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
    x       =   100,        # Number of horizontal centroids
    y       =   100,        # Number of vertical centroids
)
NV =   (
    x       =   NC.x + 1,   # Number of horizontal vertices
    y       =   NC.y + 1,   # Number of vertical vertices
)
Δ   =   (
    x   =   (abs(M.xmin)+M.xmax)/NC.x,
    y   =   (abs(M.ymin)+M.ymax)/NC.y,
)
# -------------------------------------------------------------------- #
# Erstellung des Gitters ============================================= #
x   =   (
    c       =   LinRange(M.xmin + Δ.x/2.0, M.xmax - Δ.x/2.0, NC.x),
    ce      =   LinRange(M.xmin - Δ.x/2.0, M.xmax + Δ.x/2.0, NC.x+2),
    v       =   LinRange(M.xmin, M.xmax , NV.x)    
)
y       = (
    c       =   LinRange(M.ymin + Δ.y/2.0, M.ymax - Δ.y/2.0, NC.y),
    ce      =   LinRange(M.ymin - Δ.x/2.0, M.ymax + Δ.x/2.0, NC.y+2),
    v       =   LinRange(M.ymin, M.ymax, NV.y),    
)
x1      =   ( 
    c2d     =   x.c .+ 0*y.c',
    v2d     =   x.v .+ 0*y.v', 
    vx2d    =   x.v .+ 0*y.ce',
    vy2d    =   x.ce .+ 0*y.v',
)
x   =   merge(x,x1)
y1      =   (
    c2d     =   0*x.c .+ y.c',
    v2d     =   0*x.v .+ y.v',
    vx2d    =   0*x.v .+ y.ce',
    vy2d    =   0*x.ce .+ y.v',
)
y   =   merge(y,y1)
# -------------------------------------------------------------------- #
# Tracer Advektions Verfahren ======================================== #
if FD.Method.Adv==:tracers 
    nmx,nmy     =   3,3
    noise       =   1
    nmark       =   nmx*nmy*NC.x*NC.y
    Aparam      =   :thermal
    Ma          =   IniTracer2D(Aparam,nmx,nmy,Δ,M,NC,noise)
    # RK4 weights
    rkw         =   1.0/6.0*[1.0 2.0 2.0 1.0]   # for averaging
    rkv         =   1.0/2.0*[1.0 1.0 2.0 2.0]   # for time stepping
end
# -------------------------------------------------------------------- #
# Animationssettings ================================================= #
path        =   string("./examples/AdvectionEquation/Results/")
anim        =   Plots.Animation(path, String[] )
filename    =   string("2D_advection_",Ini.T,"_",Ini.V,
                        "_",FD.Method.Adv)
save_fig    =   0
# -------------------------------------------------------------------- #
# Anfangsbedingungen ================================================= #
# Temperatur --------------------------------------------------------- #
D       =   (
    T       =   zeros(NC...),
    T_ex    =   zeros(NC.x+2,NC.y+2),
    T_exo   =   zeros(NC.x+2,NC.y+2),
    vx      =   zeros(NV.x,NV.y+1),
    vy      =   zeros(NV.x+1,NV.y),    
    vxc     =   zeros(NC...),
    vyc     =   zeros(NC...),
    vc      =   zeros(NC...),
    wt      =   zeros(NC...),
    Tmax    =   [0.0],
    Tmin    =   [0.0],
    Tmean   =   [0.0],
)
IniTemperature!(Ini.T,M,NC,Δ,D,x,y)
if FD.Method.Adv==:slf
    D.T_exo    .=  D.T_ex
end
if FD.Method.Adv==:tracers
    for k = 1:nmark
        Ma.T[k] =   FromCtoM(D.T_ex, k, Ma, x, y, Δ, NC)
    end
    CountMPC(Ma,nmark,x,y,Δ)
end
# Geschwindigkeit ---------------------------------------------------- #
IniVelocity!(Ini.V,D,NV,Δ,M,x,y)            # [ m/s ]
# Get the velocity on the centroids ---
for i = 1:NC.x, j = 1:NC.y
    D.vxc[i,j]  = (D.vx[i,j+1] + D.vx[i+1,j+1])/2
    D.vyc[i,j]  = (D.vy[i+1,j] + D.vy[i+1,j+1])/2
end
@. D.vc        = sqrt(D.vxc^2 + D.vyc^2)
# Visualize initial condition ---------------------------------------- #
if FD.Method.Adv==:tracers
    #p = scatter(Ma.x[1:Pl.Minc:end],Ma.y[1:Pl.Minc:end], 
    #        zcolor=Ma.T[1:Pl.Minc:end],color=:thermal,
    #        markersize=Pl.Msz,legend=false,colorbar=true, 
    #        aspect_ratio=:equal,xlims=(M.xmin, M.xmax), ylims=(M.ymin, M.ymax),
    #        layout=(1,2),subplot=1)
    p = heatmap(x.c,y.c,(D.T./D.Tmax)',color=:thermal, 
            aspect_ratio=:equal,xlims=(M.xmin, M.xmax), 
            ylims=(M.ymin, M.ymax),clims=(0.5, 1.0),
            colorbar=true,layout=(1,2),subplot=1)
    quiver!(p,x.c2d[1:Pl.inc:end,1:Pl.inc:end],
            y.c2d[1:Pl.inc:end,1:Pl.inc:end],
            quiver=(D.vxc[1:Pl.inc:end,1:Pl.inc:end].*Pl.sc,
                    D.vyc[1:Pl.inc:end,1:Pl.inc:end].*Pl.sc),        
            color="white",layout=(1,2),subplot=1)
    heatmap!(p,x.c,y.c,Ma.mpc',color=:inferno, 
            aspect_ratio=:equal,xlims=(M.xmin, M.xmax), ylims=(M.ymin, M.ymax),
            colorbar=true,
            layout=(1,2),subplot=2)
else
    p = heatmap(x.c , y.c, (D.T./D.Tmax)', 
            color=:thermal, colorbar=true, aspect_ratio=:equal, 
            xlabel="x", ylabel="z", 
            title="Temperature", 
            xlims=(M.xmin, M.xmax), ylims=(M.ymin, M.ymax), 
            clims=(0.5, 1.0))
    quiver!(p,x.c2d[1:Pl.inc:end,1:Pl.inc:end],y.c2d[1:Pl.inc:end,1:Pl.inc:end],
            quiver=(D.vxc[1:Pl.inc:end,1:Pl.inc:end].*Pl.sc,
                    D.vyc[1:Pl.inc:end,1:Pl.inc:end].*Pl.sc),        
            color="white")
end
if save_fig == 1
    Plots.frame(anim)
elseif save_fig == 0
    display(p)
end
# -------------------------------------------------------------------- #
# Zeit Konstanten ==================================================== #
T   =   ( 
    tmax    =   [0.0],  
    Δfac    =   1.0,    # Courant time factor, i.e. dtfac*dt_courant
    Δ       =   [0.0],
)
T.tmax[1]   =   π*((M.xmax-M.xmin)-Δ.x)/maximum(D.vc)   # t = U/v [ s ]
T.Δ[1]      =   T.Δfac * minimum((Δ.x,Δ.y)) / 
                    (sqrt(maximum(abs.(D.vx))^2 + maximum(abs.(D.vy))^2))
nt          =   ceil(Int,T.tmax[1]/T.Δ[1])
# -------------------------------------------------------------------- #
#return
# Solve advection equation ------------------------------------------- #
for i=2:nt
    display(string("Time step: ",i))    

    if FD.Method.Adv==:upwind
        upwindc2D!(D,NC,T,Δ)
    elseif FD.Method.Adv==:slf
        slfc2D!(D,NC,T,Δ)   
    elseif FD.Method.Adv==:semilag
        semilagc2D!(D,[],[],x,y,T)
    elseif FD.Method.Adv==:tracers
        # Advect tracers -------------------------------------------- #
        AdvectTracer2D(Ma,nmark,D,x,y,dt,rkw,rkv,style)
        
        # Interpolate temperature from tracers to grid -------------- #
        Markers2Cells(Ma,nmark,D.T,D.wt,x,y,Δ,Aparam)
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
    
    display(string("ΔT = ",((maximum(D.T)-D.Tmax[1])/D.Tmax[1])*100))

    # Plot Solution ---
    if mod(i,10) == 0 || i == nt
        if FD.Method.Adv==:tracers
            #p = scatter(Ma.x[1:Pl.Minc:end],Ma.y[1:Pl.Minc:end], 
            #        zcolor=Ma.T[1:Pl.Minc:end],color=:thermal,
            #        markersize=Pl.Msz,legend=false,colorbar=true, 
            #        aspect_ratio=:equal,xlims=(M.xmin, M.xmax), ylims=(M.ymin, M.ymax),
            #        layout=(1,2),subplot=1)
            p = heatmap(x.c,y.c,(D.T./D.Tmax)',color=:thermal, 
                    aspect_ratio=:equal,xlims=(M.xmin, M.xmax), 
                    ylims=(M.ymin, M.ymax),clims=(0.5, 1.0),
                    colorbar=true,layout=(1,2),subplot=1)
            quiver!(p,x.c2d[1:Pl.inc:end,1:Pl.inc:end],
                    y.c2d[1:Pl.inc:end,1:Pl.inc:end],
                    quiver=(D.vxc[1:Pl.inc:end,1:Pl.inc:end].*Pl.sc,
                            D.vyc[1:Pl.inc:end,1:Pl.inc:end].*Pl.sc),        
                    color="white",layout=(1,2),subplot=1)
            heatmap!(p,x.c,y.c,Ma.mpc',color=:inferno, 
                    aspect_ratio=:equal,xlims=(M.xmin, M.xmax), ylims=(M.ymin, M.ymax),
                    colorbar=true,
                    layout=(1,2),subplot=2)
        else
            p = heatmap(x.c , y.c, (D.T./D.Tmax)', 
                    color=:thermal, colorbar=true, aspect_ratio=:equal, 
                    xlabel="x", ylabel="z", 
                    title="Temperature", 
                    xlims=(M.xmin, M.xmax), ylims=(M.ymin, M.ymax), 
                    clims=(0.5, 1.0))
            quiver!(p,x.c2d[1:Pl.inc:end,1:Pl.inc:end],
                        y.c2d[1:Pl.inc:end,1:Pl.inc:end],
                        quiver=(D.vxc[1:Pl.inc:end,1:Pl.inc:end].*Pl.sc,
                                D.vyc[1:Pl.inc:end,1:Pl.inc:end].*Pl.sc),        
                    color="white")
        end
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
    #foreach(rm, filter(startswith(string(path,"00")), readdir(path,join=true)))
elseif save_fig == 0
    display(plot(p))
end

end # Function end

Advection_2D()
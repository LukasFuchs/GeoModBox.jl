# ----------------------------------------------------------------------- #
# Funktion zur Lösung der 1-D Advektionsgleichung unter der Annahme einer 
# konstanten horizontalen Geschwindigkeit. 
# ----------------------------------------------------------------------- #
# Vers. 1.0 - 1.11.2024
# ======================================================================= #
using Plots, Dierckx
using GeoModBox.AdvectionEquation.OneD

function Advection1D_discretization()
# Geometrische Konstanten ----------------------------------------------- #
xmin    =   0                           #   [ m ]
xmax    =   40                          #   [ m ]
# ----------------------------------------------------------------------- #
# Numerische Konstanten ------------------------------------------------- #
nc      =   100                         #   Anzahl der Gitterpunkte
Δx      =   xmax/nc                     #   Gitterlänge
ind     =   1:nc   
# ----------------------------------------------------------------------- #

xc      =   xmin+Δx/2:Δx:xmax-Δx/2      #   x-Koordinate
X       =   zeros(nc)
#xhalf   =   (xmax+abs(xmin)) / 2        

# Maximale Laufzeit des Models ------------------------------------------ #
tmax    =   40.0                        #   [ s ]
# ----------------------------------------------------------------------- #
# Horizontale Geschwindigkeit ------------------------------------------- #
vx      =   1.0                         #   [ m/s ]
# ----------------------------------------------------------------------- #
# Definition der Zeitschrittlänge --------------------------------------- #
Δtfac   =   1.0                         #   Courant-Kriterium
Δt      =   Δtfac*Δx/abs(vx)
nt      =   ceil(Int,tmax/Δt)               #   Anzahl der Zeitschritte
# ----------------------------------------------------------------------- #
# Nun können wir zum einen die FD-Methode wählen ('FTCS', 'upwind', 'downwind', 
# 'lax', 'slf', 'semi-lag' - die Gleichungen und Erläuterungen dazu sind im Detail 
# in den Folien der Vorlesung zu finden) und zum anderen das Anfangsprofil wählen 
# ('block' oder 'gaussian'):

FD          =   (Method     = (Adv=:semilag,),)
Ini         =   (T=:gaussian,)

# Animationssettings ---------------------------------------------------- #
path        =   string("./examples/AdvectionEquation/Results/")
anim        =   Plots.Animation(path, String[] )
filename    =   string("1D_advection_",Ini.T,"_",FD.Method.Adv)
save_fig    =   1
# ----------------------------------------------------------------------- #
# Tracer advection method ----------------------------------------------- #
nmx         =   3       #   Number of tracers per "cell"
# Anfangsbedingung ------------------------------------------------------ #
# Wollen wir nun die Anfangsbedingung (d.h. das Anfangstemperaturprofil) definieren: 

if Ini.T==:block 
        # Hintergrundtemperatur ----------------------------------------- #
        Tb      =   1000                #   [ K ]
        
        # Lokalität und Intensität der Temperaturanomalie --------------- #
        xTl     =   (xmax-xmin)/10
        xTr     =   xTl + (xmax-xmin)/10
        Ta      =   1500                #   [ K ]
        
        # Erstellung des Anfangstemperaturprofiles ---------------------- #
        T       =   Tb.*ones(nc)
        T[xc.>=xTl .&& xc .<= xTr ]       .=  Ta
        Tmin    =   minimum(T)
        Tmax    =   maximum(T)
        tc      =   100
elseif Ini.T==:gaussian
        # Gaußsche Temperature Verteilung ------------------------------- #
        Tb      =   1000                #   Hintergrundtempertur
        Ampl    =   500                 #   Amplitude
        sigma   =   1                   #   Standard Abweichung
        xcG     =   (xmax-xmin)/10      #   x-Koordinate des Maximums
        T       =   zeros(nc)
        @. T    =  Tb + Ampl*exp(-((xc - xcG)^2)/sigma^2)
        
        Tmin    =   minimum(T)
        Tmax    =   maximum(T)
        tc      =   Ampl
end
Told    =       zeros(nc)
Told    .=      T
Told2   =       zeros(nc)
Told2   .=      T
Tsl     =       zeros(nc)

q = plot(xc, T, xlabel="x [m]", ylabel="T [°C]", 
            title="Anfangstemperatur Verteilung", 
            markershape=:circle,label="",
            xlim=(xmin,xmax), ylim=(Tmin-10, Tmax+10))
if save_fig==0
    display(q)
end

if FD.Method.Adv==:tracers
    # Gesamtanzahl der Tracer
    nm          =   (nc)*nmx
    # Abstand der Tracer
    Δmx         =   (abs(xmin)+abs(xmax))/(nm+1)
    # x-Koordinaten der Tracer    
    xm          =   collect(xmin+Δmx:Δmx:xmax-Δmx) # .+ rand(nm).*0.5*Δmx
    # Temperatur auf den Tracern
    Tm          =   zeros(nm)
end

# Lösen der Advektionsgleichung --------------------------------------- #
for i = 2:nt
    display(string("Time step: ",i))
    
    if FD.Method.Adv==:FTCS        
        T[2:end-1]  .= 
            Told[2:end-1] .- (vx*Δt/2.0/Δx).*(Told[3:end].- T[1:end-2])
    elseif FD.Method.Adv==:upwind
        if vx > 0
            T[2:end-1] .= 
                Told[2:end-1] .- vx*Δt/Δx.*( Told[2:end-1] .- Told[1:end-2] )
        elseif vx < 0
            T[2:end-1] = 
                Told[2:end-1] .- vx*Δt/Δx.*(Told[3:end] .- T[2:end-1]) 
        end
    elseif FD.Method.Adv==:downwind
        T[2:end-1] = 
            Told[2:end-1] .- vx*Δt/Δx.*(Told[3:end] .- T[2:end-1]) 
    elseif FD.Method.Adv==:lax
        T[2:end-1] = (Told[3:end] .+ Told[1:end-2])./2 .-
            (vx*Δt/2/Δx) .* (Told[3:end] .- Told[1:end-2])
    elseif FD.Method.Adv==:slf
        if i==2
            T[2:end-1] = 
                Told[2:end-1] .- vx*Δt/Δx*(Told[3:end] .- Told[1:end-2])
        else
            T[2:end-1] = 
                Told2[2:end-1] .- vx*Δt/Δx .* (Told[3:end] .- Told[1:end-2])
        end
        Told2 .= Told
    elseif FD.Method.Adv==:semilag
        X   .=   xc .- Δt*vx
        spl     =   Spline1D(xc,T;k=1)
        Tsl     =   spl.(X)
        T       .=  Tsl      
    elseif FD.Method.Adv==:tracers 
        spl     =   Spline1D(xc,T;k=4) 
        Tm      =   spl.(xm)        
        RK4O1D!( xm, Δt, vx )
        XMT     =   hcat(xm,Tm)        
        XMT     =   sortslices(XMT,dims=1)
        xm      .=  XMT[:,1]
        Tm      .=  XMT[:,2]
        spl2    =   Spline1D(xm,Tm;k=4)
        T       =   spl2.(xc)
    end

    # Mirror boundary conditions ---------------------------------------- #
    if vx > 0
        T[end]  =    T[end-1] 
        T[1]    =    T[end]
    elseif vx < 0
        T[1]      = T[2] 
        T[end]    = T[1]
    end
    
    # Darstellung des Profils
    q = plot(xc, T, xlabel="x [m]", ylabel="T [°C]", 
            title="Anfangstemperatur Verteilung", 
            markershape=:circle,label="",
            xlim=(xmin,xmax), ylim=(Tmin-10, Tmax+10))
    if FD.Method.Adv==:tracers 
        plot!(q,xm,Tm,markershape=:pixel,label="")
    end
    if save_fig == 1
        Plots.frame(anim)
    else
        display(q)
    end
    Told        .= T
end

# Speicher Animation ---------------------------------------------------- #
if save_fig == 1
    # Write the frames to a GIF file
    Plots.gif(anim, string( path, filename, ".gif" ), fps = 15)
    foreach(rm, filter(startswith(string(path,"00")), readdir(path,join=true)))
end
# ----------------------------------------------------------------------- #
end

Advection1D_discretization()
using Plots, Dierckx, Interpolations
using GeoModBox.AdvectionEquation.OneD
function test_1D_Advection()
# Geometrische Konstanten ----------------------------------------------- #
xmin    =   0                           #   [ m ]
xmax    =   40                          #   [ m ]
# ----------------------------------------------------------------------- #
# Numerische Konstanten ------------------------------------------------- #
nc      =   100                         #   Anzahl der Gitterpunkte
Δx      =   xmax/nc                     #   Gitterlänge
ind     =   1:nc   
# ---
xc      =   xmin+Δx/2:Δx:xmax-Δx/2      #   x-Koordinate
xcwe    =   xmin-Δx/2:Δx:xmax+Δx/2      #   x-Koordinate
X       =   zeros(nc)
# ----------------------------------------------------------------------- #
# Maximale Laufzeit des Models ------------------------------------------ #
tmax    =   40.0                        #   [ s ]
# ----------------------------------------------------------------------- #
# Horizontale Geschwindigkeit ------------------------------------------- #
vx      =   1.0                         #   [ m/s ]
# ----------------------------------------------------------------------- #
# Definition der Zeitschrittlänge --------------------------------------- #
Δtfac   =   0.8                         #   Courant-Kriterium
Δt      =   Δtfac*Δx/abs(vx)
nt      =   ceil(Int,tmax/Δt)           #   Anzahl der Zeitschritte
# ----------------------------------------------------------------------- #
FD          =   (Method     = (Adv=:semilag,),)
Ini         =   (T=:gaussian,)
# Animationssettings ---------------------------------------------------- #
path        =   string("./exercise/Correction/Results/")
anim        =   Plots.Animation(path, String[] )
filename    =   string("test_",Ini.T,"_",FD.Method.Adv)
save_fig    =   0
# ----------------------------------------------------------------------- #
# Tracer advection method ----------------------------------------------- #
nmx         =   3       #   Number of tracers per "cell"
# ----------------------------------------------------------------------- #
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
TWE             =   zeros(nc+2)
TWE[2:end-1]    .=  T
TWE[1]          =   T[end]
TWE[end]        =   T[1]
TWE2            =   zeros(nc+2)

#q = plot(xc, T, xlabel="x [m]", ylabel="T [°C]", 
#        title="Anfangstemperatur Verteilung", 
#        markershape=:circle,label="",
#        xlim=(xmin,xmax), ylim=(Tmin-10, Tmax+10))
#if save_fig==0
#    display(q)
#end

if FD.Method.Adv==:tracers
    # Gesamtanzahl der Tracer
    nm          =   (nc)*nmx
    # Abstand der Tracer
    Δmx         =   (abs(xmin)+abs(xmax))/(nm+1)
    # x-Koordinaten der Tracer    
    xm          =   collect(xmin+Δmx:Δmx:xmax-Δmx) .+ rand(nm).*0.5*Δmx
    # Temperatur auf den Tracern
    Tm          =   zeros(nm)    
end
# Lösen der Advektionsgleichung --------------------------------------- #
for i = 2:nt
    display(string("Time step: ",i))
    TWE[2:end-1]    .=  T
    TWE[1]          =   T[end]
    TWE[end]        =   T[1]
    
    if FD.Method.Adv==:upwind
        if vx > 0
            T       .= 
                TWE[2:end-1] .- vx*Δt/Δx.*( TWE[2:end-1] .- TWE[1:end-2] )
        elseif vx < 0
            T .= 
                TWE[2:end-1] .- vx*Δt/Δx.*( TWE[3:end] .- TWE[2:end-1] ) 
        end
    elseif FD.Method.Adv==:semilag
        X       .=  xc .- Δt*vx
        itp_cubic   =   cubic_spline_interpolation(xcwe,TWE)
        T       .=   itp_cubic.(X)
        #spl     =   Spline1D(xcwe,TWE;k=3)
        #T       =   spl.(X)
        #spl     =   Spline1D(xcwe,TWE;k=3)
        #TWE2[2:end-1]    =   spl.(X)
        #Itp1D_Centers2Markers!(T,X,TWE,xcwe,Δx,xmin-Δx)        
    elseif FD.Method.Adv==:tracers 
        spl     =   Spline1D(xcwe,TWE;k=1) 
        Tm      =   spl.(xm)        
        #Itp1D_Centers2Markers!(Tm,xm,TWE,xcwe,Δx,xmin-Δx)
        RK4O1D!( xm, Δt, vx, xmin, xmax )
        XMT     =   hcat(xm,Tm)        
        XMT     =   sortslices(XMT,dims=1)
        xm      .=  XMT[:,1]
        Tm      .=  XMT[:,2]
        spl2    =   Spline1D(xm,Tm;k=1)
        T       =   spl2.(xc)
        #Itp1D_Markers2Centers!( T, xc, Tm, xm, Δx, xmin)
    end
    
    display(string("ΔT=",((Tmax-maximum(T))/Tmax)*100))

    # Darstellung des Profils
    q = plot(xc, T, xlabel="x [m]", ylabel="T [°C]", 
            title="Anfangstemperatur Verteilung", 
            label="")
    #plot!(q,X,T,linstyle=:dash,label="")
    if FD.Method.Adv==:tracers 
        plot!(xm,Tm,markershape=:circle,label="",linealpha=:0)        
    end    
    if save_fig == 1
        Plots.frame(anim)
    else
        display(q)
    end
end

# Speicher Animation ---------------------------------------------------- #
if save_fig == 1
    # Write the frames to a GIF file
    Plots.gif(anim, string( path, filename, ".gif" ), fps = 15)
    foreach(rm, filter(startswith(string(path,"00")), readdir(path,join=true)))
end
# ----------------------------------------------------------------------- #
end

test_1D_Advection()
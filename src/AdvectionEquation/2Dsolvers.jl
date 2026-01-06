# -------------------------------------------------------------------- #
# 2D solver for the advection equation 
# -------------------------------------------------------------------- #
using Interpolations

"""
    upwindc2D!(P,P_ex,vxc,vyc,NC,Δt,Δx,Δy)

Function to advect a property using the upwind method in two dimensions. 

    P       : 2-D array of the advected property on the centroids
    P_ex    : 2-D array including the ghost nodes
    vxc     : 2-D array of the horizontal velocity on the centroids
    vyc     : 2-D array of the vertical velocity on the centroids
    NC      : Structure or Tuple containing the number of centroids
    Δt      : Time step
    Δx      : Horizontal grid spacing
    Δy      : Vertical grid spacing
"""
function upwindc2D!(P,P_ex,vxc,vyc,NC,Δt,Δx,Δy)

    @. P     =  P_ex[2:(NC.x+1),2:(NC.y+1)] - 
            (vxc>0)*(vxc*Δt/Δx*(P_ex[2:(NC.x+1),2:(NC.y+1)] - P_ex[1:NC.x,2:(NC.y+1)])) - 
            (vxc<0)*(vxc*Δt/Δx*(P_ex[3:(NC.x+2),2:(NC.y+1)] - P_ex[2:(NC.x+1),2:(NC.y+1)])) - 
            (vyc>0)*(vyc*Δt/Δy*(P_ex[2:(NC.x+1),2:(NC.y+1)] - P_ex[2:(NC.x+1),1:NC.y])) - 
            (vyc<0)*(vyc*Δt/Δy*(P_ex[2:(NC.x+1),3:(NC.y+2)] - P_ex[2:(NC.x+1),2:(NC.y+1)]))
    # Update extende temperature field ------------------------------- #
    P_ex[2:(NC.x+1),2:(NC.y+1)]     .=  P
    # ---------------------------------------------------------------- #
end

"""
    slfc2D!

Function to advect a property using the staggered-leaped frog method in two dimensions. 

    P       : 2-D array of the advecte property on the centroids
    P_ex    : 2-D array including the ghost nodes
    P_exo   : 2-D array of the previous time step including the ghost nodes 
    vxc     : 2-D array of the horizontal velocity on the centroids
    vyc     : 2-D array of the vertical velocity on the centroids
    NC      : Structure or Tuple containing the number of centroids
    Δt      : Time step
    Δx      : Horizontal grid spacing
    Δy      : Vertical grid spacing
"""
function slfc2D!(P,P_ex,P_exo,vxc,vyc,NC,Δt,Δx,Δy)

    indx    =   2:(NC.x+1)
    indy    =   2:(NC.y+1)

    @. P  =   P_exo[indx,indy] - 
        vxc*Δt/Δx*(P_ex[indx+1,indy] - P_ex[indx.-1,indy]) - 
        vyc*Δt/Δy*(P_ex[indx,indy+1] - P_ex[indx,indy-1])
    # Update extende temperature field ------------------------------- #
    @. P_exo            =  P_ex
    P_ex[indx,indy]     .=  P
    # ---------------------------------------------------------------- #
end

"""
    semilagc2D!()

Function to advect a property using the semi-lagrange method in two dimensions. 

    P       : 2-D array of the advecte property on the centroids
    P_ex    : 2-D array including the ghost nodes
    vxc     : 2-D array of the horizontal velocity on the centroids
    vyc     : 2-D array of the vertical velocity on the centroids
    vxo     : 2-D array of the horizontal velocity of the previous time step on the centroids
    vyo     : 2-D array of the vertical velocity of the previous time step on the centroids
    x       : Structure or Tuple containing the 2-D horizontal coordinates of the centroids
    y       : Structure or Tuple containing the 2-D vertical coordinates of the centroids
    Δt      : Time step

The function calculates a mean velocity of the previous and current velocity field, if needed, 
and advects the property using a mid-point iteration scheme. 
"""
function semilagc2D!(P,P_ex,vxc,vyc,vxo,vyo,x,y,Δt)
    # mid-point iteration scheme ---
    if isempty(vxo) || isempty(vyo)
        # Im Falle das die Geschwindigkeit zeitlich konstant ist, wird die
        # aktuelle Geschwindigkeit auf die alte Geschwindigkeit
        # uebertragen.
        vxo   =   copy(vxc)
        vyo   =   copy(vyc)
    end
    
    # Mittlere Geschwindigkeit am Zentralen Punkt in der Zeit --- 
    vxcm        =   copy(vxc)
    vycm        =   copy(vyc)
    @. vxcm     =   0.5*(vxo + vxc)
    @. vycm     =   0.5*(vyo + vyc)
    
    # Initialisierung der Geschwindigkeit fuer die Iteration ---
    vxi     =   copy(vxc)
    vyi     =   copy(vyc)
    xp      =   copy(x.c2d)
    yp      =   copy(y.c2d)
    
    # Iteration ---
    for k = 1:10
        @. xp  = x.c2d - 0.5*Δt*vxi
        @. yp  = y.c2d - 0.5*Δt*vyi

        itp_linearx  =  linear_interpolation((x.c,y.c),vxcm,extrapolation_bc = Line())
        itp_lineary  =  linear_interpolation((x.c,y.c),vycm,extrapolation_bc = Line())

        vxi         .=  itp_linearx.(xp,yp)
        vyi         .=  itp_lineary.(xp,yp)
    end
    @. xp   =   x.c2d - Δt*vxi
    @. yp   =   y.c2d - Δt*vyi        

    itp_cubic   =   cubic_spline_interpolation((x.ce,y.ce),P_ex)
    P           .=  itp_cubic.(xp,yp)
        
    P_ex[2:end-1,2:end-1]   .=  P
end
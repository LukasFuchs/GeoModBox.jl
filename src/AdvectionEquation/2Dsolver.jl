# -------------------------------------------------------------------- #
# 2D solver for the advection equation 
# -------------------------------------------------------------------- #
using Interpolations

function upwindc2D!(D,NC,T,Δ)

    indx    =   2:(NC.xc+1)
    indy    =   2:(NC.yc+1)

    @. D.T     =  D.T_ex[indx,indy] - 
            (D.vxc>0)*(D.vxc*T.Δ[1]/Δ.x*(D.T_ex[indx,indy] - D.T_ex[indx-1,indy])) - 
            (D.vxc<0)*(D.vxc*T.Δ[1]/Δ.x*(D.T_ex[indx+1,indy] - D.T_ex[indx,indy])) - 
            (D.vyc>0)*(D.vyc*T.Δ[1]/Δ.y*(D.T_ex[indx,indy] - D.T_ex[indx,indy-1])) - 
            (D.vyc<0)*(D.vyc*T.Δ[1]/Δ.y*(D.T_ex[indx,indy+1] - D.T_ex[indx,indy]))
    D.T_ex[indx,indy]  .=  D.T

end

function slfc2D!(D,NC,T,Δ)

    indx    =   2:(NC.xc+1)
    indy    =   2:(NC.yc+1)

    @. D.T  =   D.T_exo[indx,indy] - 
        D.vxc*T.Δ[1]/Δ.x*(D.T_ex[indx+1,indy]-D.T_ex[indx.-1,indy]) - 
        D.vyc*T.Δ[1]/Δ.y*(D.T_ex[indx,indy+1]-D.T_ex[indx,indy-1])
    @. D.T_exo         =  D.T_ex
    D.T_ex[indx,indy]  .=  D.T
end

function semilagc2D!(D,vxo,vyo,x,y,T)
    # mid-point iteration scheme -------------------------------------------- 
    if isempty(vxo) || isempty(vyo)
        # Im Falle das die Geschwindigkeit zeitlich konstant ist, wird die
        # aktuelle Geschwindigkeit auf die alte Geschwindigkeit
        # uebertragen.
        vxo   =   copy(D.vxc)
        vyo   =   copy(D.vyc)
    end
    
    # Mittlere Geschwindigkeit am Zentralen Punkt in der Zeit --------------- 
    D.vxcm   .=   0.5.*(vxo .+ D.vxc)
    D.vycm   .=   0.5.*(vyo .+ D.vyc)
    
    # Initialisierung der Geschwindigkeit fuer die Iteration ---------------- 
    vxi     =   copy(D.vxc)
    vyi     =   copy(D.vyc)
    xp      =   copy(x.c2d)
    yp      =   copy(y.c2d)
    
    # Iteration ------------------------------------------------------------- 
    for k = 1:10
        #xp  = M.X - 0.5*dt.*vxi
        #zp  = M.Z - 0.5*dt.*vyi
        @. xp  = x.c - 0.5*T.Δ[1]*vxi
        @. yp  = y.c - 0.5*T.Δ[1]*vyi

        itp_linearx  =  linear_interpolation((x.c,y.c),D.vxc,extrapolation_bc = Line())
        itp_lineary  =  linear_interpolation((x.c,y.c),D.vyc,extrapolation_bc = Line())

        vxi         .=  itp_linearx.(xp,yp)
        vyi         .=  itp_lineary.(xp,yp)
        #vxi = interp2(M.X,M.Z,vxm,xp,zp,'linear')
        #vyi = interp2(M.X,M.Z,vzm,xp,zp,'linear')
        
        ##vxi(isnan(vxi)) = vxm(isnan(vxi));
        ##vyi(isnan(vzi)) = vzm(isnan(vzi));
    end
    @. xp   =   x.c - T.Δ[1]*vxi
    @. yp   =   y.c - T.Δ[1]*vyi
    
    itp_cubic   =   cubic_spline_interpolation((x.cew,y.cns),D.T_ex)
    D.T         .=  itp_cubic.(xp,yp)
    
    D.T_ex[2:end-1,2:end-1]  .=  D.T

    #Anew    =   interp2(M.X,M.Z,A,xp,zp,'cubic');
    
    #Anew(isnan(Anew)) = A(isnan(Anew));

    return D
end
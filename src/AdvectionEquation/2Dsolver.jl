# -------------------------------------------------------------------- #
# 2D solver for the advection equation 
# -------------------------------------------------------------------- #
function upwindc2D!(D,NC,T,Δ)

    indx    =   2:(NC.xc+1)
    indy    =   2:(NC.yc+1)

    @. D.T     =  D.T_ext[indx,indy] - 
            (D.vxc>0)*(D.vxc*T.Δ[1]/Δ.x*(D.T_ext[indx,indy] - D.T_ext[indx-1,indy])) - 
            (D.vxc<0)*(D.vxc*T.Δ[1]/Δ.x*(D.T_ext[indx+1,indy] - D.T_ext[indx,indy])) - 
            (D.vyc>0)*(D.vyc*T.Δ[1]/Δ.y*(D.T_ext[indx,indy] - D.T_ext[indx,indy-1])) - 
            (D.vyc<0)*(D.vyc*T.Δ[1]/Δ.y*(D.T_ext[indx,indy+1] - D.T_ext[indx,indy]))
    D.T_ext[indx,indy]  .=  D.T

end

function slfc2D!(D,NC,T,Δ)

    indx    =   2:(NC.xc+1)
    indy    =   2:(NC.yc+1)

    @. D.T  =   D.T_exto[indx,indy] - 
        D.vxc*T.Δ[1]/Δ.x*(D.T_ext[indx+1,indy]-D.T_ext[indx.-1,indy]) - 
        D.vyc*T.Δ[1]/Δ.y*(D.T_ext[indx,indy+1]-D.T_ext[indx,indy-1])
    @. D.T_exto         =  D.T_ext
    D.T_ext[indx,indy]  .=  D.T
end

function semilagc2D!()
    # mid-point iteration scheme -------------------------------------------- 
    if isempty(D.vxo)
        # Im Falle das die Geschwindigkeit zeitlich konstant ist, wird die
        # aktuelle Geschwindigkeit auf die alte Geschwindigkeit
        # uebertragen.
        D.vxo   =   D.vx
        D.vzo   =   D.vy
    end
    
    # Mittlere Geschwindigkeit am Zentralen Punkt in der Zeit --------------- 
    vxm     =   0.5.*(D.vxo + D.vx)
    vym     =   0.5.*(D.vzo + D.vz)
    
    # Initialisierung der Geschwindigkeit fuer die Iteration ---------------- 
    vxi     =   D.vx
    vyi     =   D.vy
    
    # Iteration ------------------------------------------------------------- 
    for k = 1:10
        xp  = M.X - 0.5*dt.*vxi
        zp  = M.Z - 0.5*dt.*vyi
        
        vxi = interp2(M.X,M.Z,vxm,xp,zp,'linear')
        vyi = interp2(M.X,M.Z,vzm,xp,zp,'linear')
        
        #vxi(isnan(vxi)) = vxm(isnan(vxi));
        #vyi(isnan(vzi)) = vzm(isnan(vzi));
    end
    xp      =   M.X - dt.*vxi;
    zp      =   M.Z - dt.*vzi;
    
    Anew    =   interp2(M.X,M.Z,A,xp,zp,'cubic');
    
    Anew(isnan(Anew)) = A(isnan(Anew));
    
end
# -------------------------------------------------------------------- #
# 2D solver for the advection equation 
# -------------------------------------------------------------------- #
function upwindc2D!(D,NC,T,Δ)

    indx    =   2:(NC.xc+1)
    indy    =   2:(NC.yc+1)

    D.T     .=  D.T_ext[indx,indy] .- 
            (D.vxc.>0).*(D.vxc.*T.Δ[1]./Δ.x.*(D.T_ext[indx,indy] .- D.T_ext[indx.-1,indy])) .- 
            (D.vxc.<0).*(D.vxc.*T.Δ[1]./Δ.x.*(D.T_ext[indx.+1,indy] .- D.T_ext[indx,indy])) .- 
            (D.vyc.>0).*(D.vyc.*T.Δ[1]./Δ.y.*(D.T_ext[indx,indy] .- D.T_ext[indx,indy.-1])) .- 
            (D.vyc.<0).*(D.vyc.*T.Δ[1]./Δ.y.*(D.T_ext[indx,indy.+1] .- D.T_ext[indx,indy]))
    D.T_ext[indx,indy]  .=  D.T

end

function slfc2D!(D,NC,T,Δ)

    indx    =   2:(NC.xc+1)
    indy    =   2:(NC.yc+1)

    D.T     .= D.T_exto[indx,indy] .- 
        D.vxc.*T.Δ[1]/Δ.x.*(D.T_ext[indx.+1,indy].-D.T_ext[indx.-1,indy]) .- 
        D.vyc.*T.Δ[1]/Δ.y.*(D.T_ext[indx,indy.+1].-D.T_ext[indx,indy.-1])
    D.T_exto            .=  D.T_ext
    D.T_ext[indx,indy]  .=  D.T
end
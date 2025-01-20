#function IniTemperature!(type,Tb,Ta,D,x)
#
#    if type==:circle
#
#    end
#
#end      
function IniVelocity!(type,D,NC,Δ,M,x,y)
    if type==:RigidBody
        # Rigid Body Rotation ---
        for i = 1:NC.xv, j = 1:NC.yv+1
            D.vx[i,j]  =    (y.cns[j]-(M.ymax-M.ymin)/2)
        end
        for i = 1:NC.xv+1, j = 1:NC.yv
            D.vy[i,j]  =   -(x.cew[i]-(M.xmax-M.xmin)/2)
        end
        
        Radx    =   zeros(size(D.vx))
        Rady    =   zeros(size(D.vy))

        @. Radx     =   sqrt((x.vx2d-(M.xmax-M.xmin)/2)^2 + (y.vx2d-(M.ymax-M.ymin)/2)^2)
        @. Rady     =   sqrt((x.vy2d-(M.xmax-M.xmin)/2)^2 + (y.vy2d-(M.ymax-M.ymin)/2)^2)

        @. D.vx[Radx>(M.xmax-M.xmin)/2-5*Δ.x]     =   0
        @. D.vy[Rady>(M.xmax-M.xmin)/2-5*Δ.x]     =   0
        
    elseif type==:ShearCell
        # Convection Cell with a Shear Deformation --- (REF?!)
        for i = 1:NC.xv, j = 1:NC.yv+1
            D.vx[i,j]   =   -sin(π*x.v[i])*cos(π*y.cns[j])
        end
        for i = 1:NC.xv+1, j = 1:NC.yv
            D.vy[i,j]   =   cos(π.*x.cew[i]).*sin(π.*y.v[j])
        end
    end
end
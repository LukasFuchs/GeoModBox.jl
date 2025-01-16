#function IniTemperature!(type,Tb,Ta,D,x)
#
#    if type==:circle
#
#    end
#
#end      
function IniVelocity!(type,D,NC,M,x,y)
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

        Radx    .=   sqrt((x.xcew2d-(xmax-xmin)/2).^2 + (Z-(zmax-zmin)/2).^2);
        #            
        #        vx(Rad>((xmax-xmin)/2)) = 0;
        #        vz(Rad>((xmax-xmin)/2)) = 0;
        #end
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
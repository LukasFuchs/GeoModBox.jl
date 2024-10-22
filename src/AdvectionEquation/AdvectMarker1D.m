function [x] = AdvectMarker1D(x,dt,vx,xmax)
% Assuming constant velocity

x1  = dt*vx; 

x2  = dt*vx; 

x3  = dt*vx; 

x4  = dt*vx; 

x   = x + (x1+2.*(x2+x3)+x4)./6; 

x(x>xmax) = x(x>xmax)-xmax; 

end
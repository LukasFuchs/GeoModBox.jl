function [XM,ZM]   = AdvectMarker2D(x,z,XM,ZM,dt,vx,vz,xmax,xmin,zmax,zmin)

% x = dt*v(n-2) --------------------------------------------------------- %
vxm = interp2(x,z,vx,XM,ZM,'linear');
vzm = interp2(x,z,vz,XM,ZM,'linear');
x1  = dt.*vxm;
z1  = dt.*vzm;

% x = dt*v(n-1)  -------------------------------------------------------- %
X   = XM + x1./2;
Z   = ZM + z1./2;
X(X>xmax) = X(X>xmax) - (max(X(X>xmax))-xmax);
X(X<xmin) = X(X<xmin) - (min(X(X<xmin))-xmin);
Z(Z>zmax) = Z(Z>zmax) - (max(Z(Z>zmax))-zmax);
Z(Z<zmin) = Z(Z<zmin) - (min(Z(Z<zmin))-zmin);
vxm = interp2(x,z,vx,X,Z,'linear');
vzm = interp2(x,z,vz,X,Z,'linear');
x2  = dt.*vxm;
z2  = dt.*vzm;

% x = dt*v(n-1)  -------------------------------------------------------- %
X   = XM + x2./2;
Z   = ZM + z2./2;
X(X>xmax) = X(X>xmax) - (max(X(X>xmax))-xmax);
X(X<xmin) = X(X<xmin) - (min(X(X<xmin))-xmin);
Z(Z>zmax) = Z(Z>zmax) - (max(Z(Z>zmax))-zmax);
Z(Z<zmin) = Z(Z<zmin) - (min(Z(Z<zmin))-zmin);
vxm = interp2(x,z,vx,X,Z,'linear');
vzm = interp2(x,z,vz,X,Z,'linear');
x3  = dt.*vxm;
z3  = dt.*vzm;

% x = dt*v(n)  ---------------------------------------------------------- %
X   = XM + x3;
Z   = ZM + z3;
X(X>xmax) = X(X>xmax) - (max(X(X>xmax))-xmax);
X(X<xmin) = X(X<xmin) - (min(X(X<xmin))-xmin);
Z(Z>zmax) = Z(Z>zmax) - (max(Z(Z>zmax))-zmax);
Z(Z<zmin) = Z(Z<zmin) - (min(Z(Z<zmin))-zmin);
vxm = interp2(x,z,vx,X,Z,'linear');
vzm = interp2(x,z,vz,X,Z,'linear');
x4  = dt.*vxm;
z4  = dt.*vzm;

XM   = XM + (x1+2.*(x2+x3)+x4)./6;
ZM   = ZM + (z1+2.*(z2+z3)+z4)./6;

XM(XM>xmax) = xmin + (XM(XM>xmax)-xmax);
XM(XM<xmin) = xmax + (XM(XM<xmin)-xmin);
ZM(ZM>zmax) = zmin + (ZM(ZM>zmax)-zmax);
ZM(ZM<zmin) = zmax + (ZM(ZM<zmin)-zmin);

% XM  = XM - min(min(XM));
% ZM  = ZM - min(min(ZM));
% nmr = size(XM,1);
%
% nz  = size(vx,1);
% nx  = size(vz,2);
%
% % dx  = max(diff(vx(1,:)));
% % dz  = max(diff(vx(:,1)));
%
% % Assuming constant velocity
% for k = 1:nmr
%     % x = dt*v(n-2) --------------------------------------------------------- %
%     [vxm,vzm] = velip(vx,vz,dx,dz,nx,nz,XM(k),ZM(k));
%     x1  = dt.*vxm;
%     z1  = dt.*vzm;
%
%     % x = dt*v(n-1)  -------------------------------------------------------- %
%     X   = XM(k) + x1./2;
%     Z   = ZM(k) + z1./2;
%     [vxm,vzm] = velip(vx,vz,dx,dz,nx,nz,X,Z);
%     x2  = dt.*vxm;
%     z2  = dt.*vzm;
%
%     % x = dt*v(n-1)  -------------------------------------------------------- %
%     X   = XM(k) + x2./2;
%     Z   = ZM(k) + z2./2;
%     [vxm,vzm] = velip(vx,vz,dx,dz,nx,nz,X,Z);
%     x3  = dt.*vxm;
%     z3  = dt.*vzm;
%
%     % x = dt*v(n)  ---------------------------------------------------------- %
%     X   = XM(k) + x3;
%     Z   = ZM(k) + z3;
%     [vxm,vzm] = velip(vx,vz,dx,dz,nx,nz,X,Z);
%     x4  = dt.*vxm;
%     z4  = dt.*vzm;
%
%     XM(k)   = XM(k) + (x1+2.*(x2+x3)+x4)./6;
%     ZM(k)   = ZM(k) + (z1+2.*(z2+z3)+z4)./6;
%
% end

% end


end
%
% function [uip,wip] = velip(u,w,dx,dz,nx,nz,x,z)
% %% VELIP
% %
% % ----------------------------------------------------------------------- %
%
% i   = round(x/dx) + 1;
% j   = round(z/dz) + 1;
%
% % x,z outside the box? Then extrapolation from closest grid cell
% if(i>=nx)
%     i = nx-1;
% end
% if(i<1)
%     i = 1;
% end
% if(j>=nz)
%     j = nz-1;
% end
% if(j<1)
%     j = 1;
% end
%
% i1      = i+1;
% j1      = j+1;
%
% cx      = (x-dx*(i-1))/dx;
% cz      = (z-dz*(j-1))/dz;
%
% if(i1 > nx)
%     i1 = nx;
% end
% if(j1 > nz)
%     j1 = nz;
% end
%
% u1      = u(j,i) + cx*(u(j,i1) - u(j,i));
% u2      = u(j,i1) + cx*(u(j1,i1) - u(j1,i));
% uip     = u1 + cz*(u2-u1);
%
% w1      = w(j,i) + cx*(w(j,i1) - w(j,i));
% w2      = w(j,i1) + cx*(w(j1,i1) - w(j1,i));
% wip     = w1 + cz*(w2-w1);
%
% % end
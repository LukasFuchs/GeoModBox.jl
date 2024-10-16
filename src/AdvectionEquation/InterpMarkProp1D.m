function [] = InterpMarkProp1D(Bg,Bm,Xm,Xg)

nx  = length(Xg);       % Number of grid points

dx  = min(diff(Xg)); 

for i = 1:nx
    dxm     = x(i) - Xm; 
    
    wm      = 1 - dxm/dx;
    
    Bg(i)   = sum(Bm.*wm)/sum(wm); 
end

end
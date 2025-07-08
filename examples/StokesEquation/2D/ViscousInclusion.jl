using GeoModBox, Printf, Plots
using GeoModBox.InitialCondition, GeoModBox.Tracers.TwoD
using GeoModBox.MomentumEquation.TwoD
using Base.Threads, LinearAlgebra

function ViscousInclusion()
## ===================== Some initial definitions ======================= #
save_fig    =   0
plotfields  =:yes
#Py.scale        =   'yes';
chkvel      =   1; 
FDMethod    =:dc
# ----------------------------------------------------------------------- #
# Plot Settings ========================================================= #
Pl  =   (
    qinc    =   10,
    mainc   =   2,
    qsc     =   5e-4# 100*(60*60*24*365.25)
)
# ------------------------------------------------------------------- #
# Define Initial Condition ============================================== #
Ini         =   (
    V=:SimpleShear,
    p=:Inclusion,
) 
# inclusions bedingungen
α           =   0;             # positive -> counter clockwise
EllA        =   2e2 # 1.75e3; [m]
EllB        =   2e2 # 0.25e3; [m] 
# ----------------------------------------------------------------------- #
## ==================== Define model geometry constants ================= #
M       =   Geometry(
        ymin    =   -1.0e3,     #   Model depth [ m ]
        ymax    =   0.0,        
        xmin    =   0.0,
        xmax    =   1.0e3,      #   Model length [ m ]
)
# ----------------------------------------------------------------------- #
## ====================== Define the numerical grid ===================== #
NC  =   ( 
    x   =   150, 
    y   =   150, 
)
NV  =   (
    x   =   NC.x + 1,
    y   =   NC.y + 1,
)
Δ       =   GridSpacing(
    x   =   (M.xmax - M.xmin)/NC.x,
    y   =   (M.ymax - M.ymin)/NC.y,
)
x       =   (
    c   =   LinRange(M.xmin+Δ.x/2,M.xmax-Δ.x/2,NC.x),
    ce  =   LinRange(M.xmin - Δ.x/2.0, M.xmax + Δ.x/2.0, NC.x+2),
    v   =   LinRange(M.xmin,M.xmax,NV.x),
)
y       =   (
    c   =   LinRange(M.ymin+Δ.y/2,M.ymax-Δ.y/2,NC.y),
    ce  =   LinRange(M.ymin - Δ.x/2.0, M.ymax + Δ.x/2.0, NC.y+2),
    v   =   LinRange(M.ymin,M.ymax,NV.y),
)
x1      =   (
    c2d     =   x.c .+ 0*y.c',
    v2d     =   x.v .+ 0*y.v', 
    vx2d    =   x.v .+ 0*y.ce',
    vy2d    =   x.ce .+ 0*y.v',
)
x   =   merge(x,x1)
y1      =   (
    c2d     =   0*x.c .+ y.c',
    v2d     =   0*x.v .+ y.v',
    vx2d    =   0*x.v .+ y.ce',
    vy2d    =   0*x.ce .+ y.v',
)
y   =   merge(y,y1)
# ----------------------------------------------------------------------- #
## ====================== Define physical constants ===================== #
g       =   10.0                #   Gravitational acceleration [ m/s^2 ]
# 0 - upper layer; 1 - lower layer
η₀      =   1e19                #   Viscosity composition 0 [ Pa s ]
η₁      =   1e23                #   Viscosity composition 1 - Inclusion [ Pa s]
# ηᵣ      =   log10(η₁/η₀)
η       =   [η₀,η₁]             #   Viscosity for phases 

ρ₀      =   3200.0              #   Density composition 0 [ kg/m^3 ]
ρ₁      =   3200.0              #   Density composition 1 [ kg/m^3 ]
ρ       =   [ρ₀,ρ₁]             #   Density for phases

phase   =   [0,1]
# ----------------------------------------------------------------------- #
# Allocation ======================================================== #
D       =   (
    ρ       =   zeros(Float64,(NC...)),
    p       =   zeros(Float64,(NC...)),
    cp      =   zeros(Float64,(NC...)),
    vx      =   zeros(Float64,(NV.x,NV.y+1)),
    vy      =   zeros(Float64,(NV.x+1,NV.y)),    
    Pt      =   zeros(Float64,(NC...)),
    vxc     =   zeros(Float64,(NC...)),
    vyc     =   zeros(Float64,(NC...)),
    vc      =   zeros(Float64,(NC...)),
    wt      =   zeros(Float64,(NC.x,NC.y)),
    wtv     =   zeros(Float64,(NV.x,NV.y)),
    ηc      =   zeros(Float64,NC...),
    ηv      =   zeros(Float64,NV...),
)
# ------------------------------------------------------------------- #
# Needed for the defect correction solution ---
divV        =   zeros(Float64,NC...)
ε           =   (
    xx      =   zeros(Float64,NC...), 
    yy      =   zeros(Float64,NC...), 
    xy      =   zeros(Float64,NV...),
)
τ           =   (
    xx      =   zeros(Float64,NC...), 
    yy      =   zeros(Float64,NC...), 
    xy      =   zeros(Float64,NV...),
)
# ----------------------------------------------------------------------- #
# Boundary Conditions =================================================== #
VBC     =   (
    type    =   (E=:const,W=:const,S=:const,N=:const),
    val     =   (E=zeros(NV.y),W=zeros(NV.y),S=zeros(NV.x),N=zeros(NV.x),
                vxE=zeros(NC.y),vxW=zeros(NC.y),vyS=zeros(NC.x),vyN=zeros(NC.x)),
)
# ----------------------------------------------------------------------- #
# Initial Condition ===================================================== #
IniVelocity!(Ini.V,D,VBC,NC,NV,Δ,M,x,y;rad=EllA,mus_i=η₁/η₀)
for i = 1:NC.x
    for j = 1:NC.y
        D.vxc[i,j]  = (D.vx[i,j+1] + D.vx[i+1,j+1])/2
        D.vyc[i,j]  = (D.vy[i+1,j] + D.vy[i+1,j+1])/2
    end
end
@. D.vc        = sqrt(D.vxc^2 + D.vyc^2)
@show minimum(D.vx),maximum(D.vx), minimum(D.vy),maximum(D.vy)
vxana       =   copy(D.vx)
vyana       =   copy(D.vy)
Pta         =   copy(D.Pt)
# ----------------------------------------------------------------------- #
# Tracer Advection ====================================================== #
nmx,nmy     =   5,5
noise       =   0
nmark       =   nmx*nmy*NC.x*NC.y
Aparam      =   :phase
MPC         =   (
    c       =   zeros(Float64,(NC.x,NC.y)),
    v       =   zeros(Float64,(NV.x,NV.y)),
    th      =   zeros(Float64,(nthreads(),NC.x,NC.y)),
    thv     =   zeros(Float64,(nthreads(),NV.x,NV.y)),
)
MPC1        = (
    PG_th   =   [similar(D.ρ) for _ = 1:nthreads()],    # per thread
    PV_th   =   [similar(D.ηv) for _ = 1:nthreads()],   # per thread
    wt_th   =   [similar(D.wt) for _ = 1:nthreads()],   # per thread
    wtv_th  =   [similar(D.wtv) for _ = 1:nthreads()],  # per thread
)
MPC     =   merge(MPC,MPC1)
Ma      =   IniTracer2D(Aparam,nmx,nmy,Δ,M,NC,noise,Ini.p,phase;
                ellA=EllA,ellB=EllB,α=α)
# RK4 weights ---
rkw     =   1.0/6.0*[1.0 2.0 2.0 1.0]   # for averaging
rkv     =   1.0/2.0*[1.0 1.0 2.0 2.0]   # for time stepping
# Count marker per cell ---
CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,NV,1)
# Interpolate from markers to cell ---
Markers2Cells(Ma,nmark,MPC.PG_th,D.ρ,MPC.wt_th,D.wt,x,y,Δ,Aparam,ρ)
Markers2Cells(Ma,nmark,MPC.PG_th,D.p,MPC.wt_th,D.wt,x,y,Δ,Aparam,phase)
Markers2Vertices(Ma,nmark,MPC.PV_th,D.ηv,MPC.wtv_th,D.wtv,x,y,Δ,Aparam,η)
@. D.ηc     =   0.25 * (D.ηv[1:end-1,1:end-1] + 
                        D.ηv[2:end-0,1:end-1] + 
                        D.ηv[1:end-1,2:end-0] + 
                        D.ηv[2:end-0,2:end-0])
# ------------------------------------------------------------------- #
# System of Equations =============================================== #
# Iterations
niter   =   10
ϵ       =   1e-8
# Numbering, without ghost nodes! ---
off    = [  NV.x*NC.y,                          # vx
            NV.x*NC.y + NC.x*NV.y,              # vy
            NV.x*NC.y + NC.x*NV.y + NC.x*NC.y]  # Pt

Num    =    (
    Vx  =   reshape(1:NV.x*NC.y, NV.x, NC.y), 
    Vy  =   reshape(off[1]+1:off[1]+NC.x*NV.y, NC.x, NV.y), 
    Pt  =   reshape(off[2]+1:off[2]+NC.x*NC.y,NC...),
            )
δx      =   zeros(maximum(Num.Pt))
F       =   zeros(maximum(Num.Pt))
χ       =   zeros(maximum(Num.Pt))      #   Unknown Vector ME
rhs     =   zeros(maximum(Num.Pt))      #   Right-hand Side M
# Residuals ---
Fm     =    (
    x       =   zeros(Float64,NV.x, NC.y), 
    y       =   zeros(Float64,NC.x, NV.y)
)
FPt     =   zeros(Float64,NC...)      
# ----------------------------------------------------------------------- #
# ## ========================== Scale Parameters ========================== #
# switch lower(Py.scale)
#     case 'yes'
#         [M,N,D,T,S]         =   ScaleParameters(B,M,Py,N,D,T);
# end
# # ----------------------------------------------------------------------- #
# Solution ============================================================== #
# Momentum Equation ===
# Initial Residual ---------------------------------------------- #
D.vx[2:end-1,:]     .=  0.0
D.vy[:,2:end-1]     .=  0.0
D.Pt    .=  1.0
if FDMethod==:dc
    for iter=1:niter
        Residuals2D!(D,VBC,ε,τ,divV,Δ,D.ηc,D.ηv,g,Fm,FPt)
        F[Num.Vx]   =   Fm.x[:]
        F[Num.Vy]   =   Fm.y[:]
        F[Num.Pt]   =   FPt[:]
        @printf("||R|| = %1.4e\n", norm(F)/length(F))
        norm(F)/length(F) < ϵ ? break : nothing
        # Assemble Coefficients ========================================= #
        K       =   Assembly(NC, NV, Δ, D.ηc, D.ηv, VBC, Num)
        # --------------------------------------------------------------- #
        # Solution of the linear system ================================= #
        δx      =   - K \ F
        # --------------------------------------------------------------- #
        # Update Unknown Variables ====================================== #
        D.vx[:,2:end-1]     .+=  δx[Num.Vx]
        D.vy[2:end-1,:]     .+=  δx[Num.Vy]
        D.Pt                .+=  δx[Num.Pt]
    end
elseif FDMethod==:direct
    # Update K ---
    K   =   Assembly( NC, NV, Δ, D.ηc, D.ηv, VBC, Num )
    # Update RHS ---
    # rhs term defined by the Boussinesq approximation
    rhs     =   updaterhs( NC, NV, Δ, D.ηc, D.ηv, D.ρ, -g, VBC, Num )
    # Solve System of Equations ---
    χ       =   K \ rhs
    # Update Unknown Variables ---
    D.vx[:,2:end-1]     .=  χ[Num.Vx]
    D.vy[2:end-1,:]     .=  χ[Num.Vy]
    D.Pt                .=  χ[Num.Pt]
end
# --------------------------------------------------------------- #
# Get the velocity on the centroids ---
for i = 1:NC.x
    for j = 1:NC.y
        D.vxc[i,j]  = (D.vx[i,j+1] + D.vx[i+1,j+1])/2
        D.vyc[i,j]  = (D.vy[i+1,j] + D.vy[i+1,j+1])/2
    end
end
@. D.vc        = sqrt(D.vxc^2 + D.vyc^2)

p = heatmap(x.v./1e3,y.ce./1e3,D.vx',color=:berlin,
            xlabel="x[km]",ylabel="y[km]",colorbar=true,
            title="v_x",
            aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),                             
            ylims=(M.ymin/1e3, M.ymax/1e3),
            layout=(2,2),subplot=1)
heatmap!(p,x.ce./1e3,y.v./1e3,D.vy',color=:berlin,
            xlabel="x[km]",ylabel="y[km]",colorbar=true,
            title="v_y",
            aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),                             
            ylims=(M.ymin/1e3, M.ymax/1e3),
            layout=(2,2),subplot=2)
heatmap!(p,x.c./1e3,y.c./1e3,D.Pt',color=:glasgow,
            xlabel="x[km]",ylabel="y[km]",colorbar=true,
            title="P_t",
            aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),                             
            ylims=(M.ymin/1e3, M.ymax/1e3),
            layout=(2,2),subplot=3)
heatmap!(p,x.c./1e3,y.c./1e3,log10.(D.ηc'),color=:lipari10,
            xlabel="x[km]",ylabel="y[km]",colorbar=true,
            title="phase",
            aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),                             
            ylims=(M.ymin/1e3, M.ymax/1e3),
            layout=(2,2),subplot=4)
quiver!(p,x.c2d[1:Pl.qinc:end,1:Pl.qinc:end]./1e3,
            y.c2d[1:Pl.qinc:end,1:Pl.qinc:end]./1e3,
            quiver=(D.vxc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc,
                    D.vyc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc),        
                    la=0.5,color="white",layout=(2,2),subplot=4)
# p = heatmap(x.v./1e3,y.ce./1e3,D.vx',color=:berlin,
#             xlabel="x[km]",ylabel="y[km]",colorbar=true,
#             title="v_x",
#             aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),                             
#             ylims=(M.ymin/1e3, M.ymax/1e3),
#             layout=(1,3),subplot=1)
# heatmap!(p,x.ce./1e3,y.v./1e3,D.vy',color=:berlin,
#             xlabel="x[km]",ylabel="y[km]",colorbar=true,
#             title="v_y",
#             aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),                             
#             ylims=(M.ymin/1e3, M.ymax/1e3),
#             layout=(1,3),subplot=2)
# heatmap!(p,x.c./1e3,y.c./1e3,D.Pt',color=:glasgow,
#             xlabel="x[km]",ylabel="y[km]",colorbar=true,
#             title="|v|",
#             aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),                             
#             ylims=(M.ymin/1e3, M.ymax/1e3),
#             layout=(1,3),subplot=3)
# heatmap!(x.v./1e3,y.ce./1e3,vxana',color=:berlin,
#             xlabel="x[km]",ylabel="y[km]",colorbar=true,
#             title="v_x_ana",
#             aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),                             
#             ylims=(M.ymin/1e3, M.ymax/1e3),
#             layout=(3,3),subplot=2)
# heatmap!(p,x.ce./1e3,y.v./1e3,vyana',color=:berlin,
#             xlabel="x[km]",ylabel="y[km]",colorbar=true,
#             title="v_y",
#             aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),                             
#             ylims=(M.ymin/1e3, M.ymax/1e3),
#             layout=(3,3),subplot=5)
# heatmap!(p,x.c./1e3,y.c./1e3,Pta',color=:glasgow,
#             xlabel="x[km]",ylabel="y[km]",colorbar=true,
#             title="|v|",
#             aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),                             
#             ylims=(M.ymin/1e3, M.ymax/1e3),
#             layout=(3,3),subplot=8)
# quiver!(p,x.c2d[1:Pl.qinc:end,1:Pl.qinc:end]./1e3,
#             y.c2d[1:Pl.qinc:end,1:Pl.qinc:end]./1e3,
#             quiver=(D.vxc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc,
#                     D.vyc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc),        
#                     la=0.5,color="white",layout=(2,2),subplot=3)
# heatmap!(p,x.c./1e3,y.c./1e3,D.p',color=:lipari10,
#             xlabel="x[km]",ylabel="y[km]",colorbar=true,
#             title="phase",
#             aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),                             
#             ylims=(M.ymin/1e3, M.ymax/1e3),
#             layout=(2,2),subplot=4)

display(p)
# ## ========================= Time loop ================================= ##
#     ## =========== Interpolate velocity onto the regular grid =========== #
#     [ID]        =   InterpStaggered(D,ID,N,'velocity');
#     D.meanV(it) =   rms(ID.vx(:) + ID.vz(:));
#     [ID]        =   GetStrainRate(ID,N);
#     ID.tauII    =   ID.eII.*D.eta.*2;
#     ID.psi      =   ID.eII.*ID.tauII;
#     # =================================================================== #
#     [D]         =   GetStrainRateStag(D,N);
#     [D]         =   GetStressStag(D,N);
#     D.psi       =   D.eII.*D.tauII;
        
#     incind      =   D.C > 1.5; 
#     #             matind      =   D.C <= 1.5;   
#     # =================================================================== #
#     ## ========================== Plot data ============================= #
#     Pl.time     =   '';
#     Pl.xlab     =   '$$x$$';
#     Pl.zlab     =   '$$z$$';
    
#     if (mod(it,5)==0||it==1)
#         figure(1) # ----------------------------------------------------- #
#         clf
#         ax1 = subplot(2,2,1);
#         D.eta(~incind) = NaN;
#         plotfield(log10(D.eta),M.X,M.Z,Pl,'pcolor',...
#             '$$log_{10} (\ \eta\ )$$','quiver',ID.vx,ID.vz);
#         colormap(ax1,flipud(Pl.lapaz))
#         ax2 = subplot(2,2,2);
#         D.psi(~incind) = NaN;
#         plotfield((D.psi),M.X,M.Z,Pl,'pcolor',...
#             '$$log_{10} (\ \psi\ )$$');
#         colormap(ax2,Pl.imola)
#         ax3 = subplot(2,2,3);
#         D.eII(~incind) = NaN;
#         plotfield((D.eII),M.X,M.Z,Pl,'pcolor',...
#             '$$log_{10} (\ \dot\varepsilon_{II}\ )$$')
#         colormap(ax3,Pl.batlowW)
#         ax4 = subplot(2,2,4);
#         D.tauII(~incind) = NaN;
#         plotfield((D.tauII),M.X,M.Z,Pl,'pcolor',...
#             '$$log_{10} (\ \tau_{II}\ )$$')
#         colormap(ax4,Pl.nuuk)
# #         figure(3) # ----------------------------------------------------- #
# #         clf
# #         ax1 = subplot(2,2,1);
# #         D.eta(~incind) = NaN;
# #         plotfield(log10(D.eta),M.X,M.Z,Pl,'pcolor',...
# #             '$$log_{10} (\ \eta\ )$$','quiver',ID.vx,ID.vz);
# #         colormap(ax1,flipud(Pl.lapaz))
# #         ax2 = subplot(2,2,2);
# #         ID.psi(~incind) = NaN;
# #         plotfield((ID.psi),M.X,M.Z,Pl,'pcolor',...
# #             '$$log_{10} (\ \psi\ )$$');
# #         colormap(ax2,Pl.imola)
# #         ax3 = subplot(2,2,3);
# #         ID.eII(~incind) = NaN;
# #         plotfield((ID.eII),M.X,M.Z,Pl,'pcolor',...
# #             '$$log_{10} (\ \dot\varepsilon_{II}\ )$$')
# #         colormap(ax3,Pl.batlowW)
# #         ax4 = subplot(2,2,4);
# #         ID.tauII(~incind) = NaN;
# #         plotfield((ID.tauII),M.X,M.Z,Pl,'pcolor',...
# #             '$$log_{10} (\ \tau_{II}\ )$$')
# #         colormap(ax4,Pl.nuuk)
#     end       
#     # =================================================================== #
#     if B.chkvel == 1
#         ID  =   CheckContinuum(ID,N,M,Ma,Pl);
#     end
    
# end
# # Staggered grid coordinates
# [M.xVx,M.zVx] = meshgrid(M.x,M.z1);
# [M.xVz,M.zVz] = meshgrid(M.x1,M.z);

# D.Pe        = abs(D.P(2:end,2:end)-D.Pa);
# D.vxe       = abs(D.vx(1:end-1,:)-D.Vxa);
# D.vze       = abs(D.vz(:,1:end-1)-D.Vza);

# figure(2)
# ax1 = subplot(3,3,1); plotfield(D.Vxa,M.xVx,M.zVx,Pl,...
#     'pcolor','$$vx\ (analytical)$$')
# colormap(ax1,Pl.imola)
# ax2 = subplot(3,3,2); plotfield(D.vx(1:end-1,:),M.xVx,M.zVx,Pl,...
#     'pcolor','$$vx\ (numerical)$$')
# colormap(ax2,Pl.imola)
# ax3 = subplot(3,3,3); plotfield(D.vxe,M.xVx,M.zVx,Pl,...
#     'pcolor','$$vx\ (error)$$')
# colormap(ax3,Pl.batlowW)
# ax4 = subplot(3,3,4); plotfield(D.Vza,M.xVz,M.zVz,Pl,...
#     'pcolor','$$vz\ (analytical)$$')
# colormap(ax4,Pl.imola)
# ax5 = subplot(3,3,5); plotfield(D.vz(:,1:end-1),M.xVz,M.zVz,Pl,...
#     'pcolor','$$vz\ (numerical)$$')
# colormap(ax5,Pl.imola)
# ax6 = subplot(3,3,6); plotfield(D.vze,M.xVz,M.zVz,Pl,...
#     'pcolor','$$vz\ (error)$$')
# colormap(ax6,Pl.batlowW)
# ax7 = subplot(3,3,7); plotfield(D.Pa,M.X1,M.Z1,Pl,...
#     'pcolor','$$P\ (analytical)$$')
# colormap(ax7,Pl.hawaii)     
# ax8 = subplot(3,3,8); plotfield(D.P(2:end,2:end),M.X1,M.Z1,Pl,...
#     'pcolor','$$P\ (numerical)$$')
# colormap(ax8,Pl.hawaii)
# ax9 = subplot(3,3,9); plotfield(D.Pe,M.X1,M.Z1,Pl,...
#     'pcolor','$$P\ (error)$$')
# colormap(ax9,Pl.batlowW)

# psiinc1 = mean(ID.psi(incind));
# psiinc2 = -sum(ID.psi(incind))./sum(incind(:).*N.dx.*N.dz);
# psiinc3 = sum(ID.psi(incind))./pi/B.EllA/B.EllB;

# ----------------------------------------------------------------------- #
# =============================== END =================================== #
# ----------------------------------------------------------------------- #
end

ViscousInclusion()
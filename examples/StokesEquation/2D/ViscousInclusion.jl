using GeoModBox, Printf, Plots
using GeoModBox.InitialCondition, GeoModBox.Tracers.TwoD
using GeoModBox.MomentumEquation.TwoD
using Base.Threads, LinearAlgebra

function ViscousInclusion()
## ===================== Some initial definitions ======================= #
save_fig    =   0
plotfields  =:yes
#Py.scale        =   'yes';
FDMethod    =:dc
# ----------------------------------------------------------------------- #
# Plot Settings ========================================================= #
Pl  =   (
    qinc    =   8,
    mainc   =   2,
    qsc     =   1e8# 100*(60*60*24*365.25)
)
# ------------------------------------------------------------------- #
# Define Initial Condition ============================================== #
Ini         =   (
    V=:PureShear,
    p=:Inclusion,
    ε = 1e-12,
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
    x   =   100, 
    y   =   100, 
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
IniVelocity!(Ini.V,D,VBC,NC,NV,Δ,M,x,y;Ini.ε)
# for i = 1:NC.x
#     for j = 1:NC.y
#         D.vxc[i,j]  = (D.vx[i,j+1] + D.vx[i+1,j+1])/2
#         D.vyc[i,j]  = (D.vy[i+1,j] + D.vy[i+1,j+1])/2
#     end
# end
# @. D.vc        = sqrt(D.vxc^2 + D.vyc^2)
# @show minimum(D.vx),maximum(D.vx), minimum(D.vy),maximum(D.vy)
# vxana       =   copy(D.vx)
# vyana       =   copy(D.vy)
# Pta         =   copy(D.Pt)
# Get analytical Solution
AnaSol  =   (
    Pa      =   zeros(Float64,NC...),
    Vxa     =   zeros(Float64,NV.x,NC.y),
    Vya     =   zeros(Float64,NC.x,NV.y),
    Vx_N    =   zeros(Float64,NV.x,1),
    Vx_S    =   zeros(Float64,NV.x,1),
    Vx_W    =   zeros(Float64,NC.y,1),
    Vx_E    =   zeros(Float64,NC.y,1),
    Vy_N    =   zeros(Float64,NC.x,1),
    Vy_S    =   zeros(Float64,NC.x,1),
    Vy_W    =   zeros(Float64,NV.y,1),
    Vy_E    =   zeros(Float64,NV.y,1),
)

Dani_Solution_vec!(Ini.V,AnaSol,M,x,y,EllA,η₁/η₀,NC,NV)

# Boundary Conditions ---
# Horizontal velocity 
VBC.val.S    .=  AnaSol.Vx_S./1e12
VBC.val.N    .=  AnaSol.Vx_N./1e12
VBC.val.vxE  .=  AnaSol.Vx_E./1e12
VBC.val.vxW  .=  AnaSol.Vx_W./1e12

# Vertical velocity 
VBC.val.E    .=  AnaSol.Vy_E./1e12
VBC.val.W    .=  AnaSol.Vy_W./1e12
VBC.val.vyS  .=  AnaSol.Vy_S./1e12
VBC.val.vyN  .=  AnaSol.Vy_N./1e12

@. D.vx[1,2:end-1]      =   AnaSol.Vx_W/1e12
@. D.vx[end,2:end-1]    =   AnaSol.Vx_E/1e12
@. D.vy[2:end-1,1]      =   AnaSol.Vy_S/1e12
@. D.vy[2:end-1,end]    =   AnaSol.Vy_N/1e12
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
# # Count marker per cell ---
# CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,NV,1)
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
niter   =   50
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
@. D.vx     *=  1e12
@. D.vy     *=  1e12
@. D.Pt     *=  1e12
# # Get the velocity on the centroids ---
# for i = 1:NC.x
#     for j = 1:NC.y
#         D.vxc[i,j]  = (D.vx[i,j+1] + D.vx[i+1,j+1])/2
#         D.vyc[i,j]  = (D.vy[i+1,j] + D.vy[i+1,j+1])/2
#     end
# end
# @. D.vc        = sqrt(D.vxc^2 + D.vyc^2)
Pe      =   copy(D.Pt)
@. Pe   =   abs((D.Pt-AnaSol.Pa))
Vxe     =   copy(AnaSol.Vxa)
@. Vxe  =   abs((D.vx[:,2:end-1]-AnaSol.Vxa))
Vye     =   copy(AnaSol.Vya)
@. Vye  =   abs((D.vy[2:end-1,:]-AnaSol.Vya))

p = heatmap(x.v./1e3,y.ce./1e3,(D.vx)',color=:berlin,
            xlabel="x[km]",ylabel="y[km]",colorbar=true,
            title="Numerical",
            aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),                             
            ylims=(M.ymin/1e3, M.ymax/1e3),
            layout=(3,3),subplot=1)
heatmap!(p,x.ce./1e3,y.v./1e3,(D.vy)',color=:berlin,
            xlabel="x[km]",ylabel="y[km]",colorbar=true,
            title="v_y [cm/a]",
            aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),                             
            ylims=(M.ymin/1e3, M.ymax/1e3),
            layout=(3,3),subplot=4)
heatmap!(p,x.c./1e3,y.c./1e3,D.Pt',color=:glasgow,
            xlabel="x[km]",ylabel="y[km]",colorbar=true,
            title="P_t",
            aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),                             
            ylims=(M.ymin/1e3, M.ymax/1e3),
            layout=(3,3),subplot=7)
heatmap!(p,x.v./1e3,y.c./1e3,(AnaSol.Vxa)',color=:berlin,
            xlabel="x[km]",ylabel="y[km]",colorbar=true,
            title="Analytical",
            aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),                             
            ylims=(M.ymin/1e3, M.ymax/1e3),
            layout=(3,3),subplot=2)
heatmap!(p,x.c./1e3,y.v./1e3,(AnaSol.Vya)',color=:berlin,
            xlabel="x[km]",ylabel="y[km]",colorbar=true,
            title="v_y [cm/a]",
            aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),                             
            ylims=(M.ymin/1e3, M.ymax/1e3),
            layout=(3,3),subplot=5)
heatmap!(p,x.c./1e3,y.c./1e3,AnaSol.Pa',color=:glasgow,
            xlabel="x[km]",ylabel="y[km]",colorbar=true,
            title="P_t",
            aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),                             
            ylims=(M.ymin/1e3, M.ymax/1e3),
            layout=(3,3),subplot=8)
heatmap!(p,x.v./1e3,y.c./1e3,(Vxe)',color=:berlin,
            xlabel="x[km]",ylabel="y[km]",colorbar=true,
            title="Analytical",
            aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),                             
            ylims=(M.ymin/1e3, M.ymax/1e3),
            layout=(3,3),subplot=3)
heatmap!(p,x.c./1e3,y.v./1e3,(Vye)',color=:berlin,
            xlabel="x[km]",ylabel="y[km]",colorbar=true,
            title="v_y [cm/a]",
            aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),                             
            ylims=(M.ymin/1e3, M.ymax/1e3),
            layout=(3,3),subplot=6)
heatmap!(p,x.c./1e3,y.c./1e3,Pe',color=:glasgow,
            xlabel="x[km]",ylabel="y[km]",colorbar=true,
            title="P_t",
            aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),                             
            ylims=(M.ymin/1e3, M.ymax/1e3),
            layout=(3,3),subplot=9)


# heatmap!(p,x.c./1e3,y.c./1e3,log10.(D.ηc'),color=reverse(cgrad(:roma)),
#             xlabel="x[km]",ylabel="y[km]",colorbar=true,
#             title="η",
#             aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),                             
#             ylims=(M.ymin/1e3, M.ymax/1e3),
#             layout=(2,2),subplot=4)
# quiver!(p,x.c2d[1:Pl.qinc:end,1:Pl.qinc:end]./1e3,
#             y.c2d[1:Pl.qinc:end,1:Pl.qinc:end]./1e3,
#             quiver=(D.vxc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc,
#                     D.vyc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc),        
#                     la=0.5,color="white",layout=(2,2),subplot=4)

display(p)

# ----------------------------------------------------------------------- #
# =============================== END =================================== #
# ----------------------------------------------------------------------- #
end

ViscousInclusion()
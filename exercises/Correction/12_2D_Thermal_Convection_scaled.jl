using Plots, ExtendableSparse
using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
using GeoModBox.AdvectionEquation.TwoD, GeoModBox.HeatEquation.TwoD
# using GeoModBox.Tracers.TwoD
# using Base.Threads
using Statistics, LinearAlgebra
using Printf

function ThermalConvection_scaled()
# Define numerical methods ========================================== #
# Advection Scheme ---
#   1) upwind, 2) slf, 3) semilag, 4) tracers
#       --- slf instable , tracers need to be modified ---
# Diffusion Scheme --- 
#   1) explicit, 2) implicit, 3) CNA, 4) ADI, 5) dc
#       dc - source term missing!
# Momentum Equation --- 
#   1) direct, 2) dc 
FD          =   (Method     = (
    Diff=:explicit,
    Adv=:upwind,
    Mom=:direct),
)
# Define Initial Condition ---
# Temperature - 
#   1) circle, 2) gaussian, 3) block, 4) linear, 5) lineara
# !!! Gaussian is not working!!! 
# Velocity - 
#   1) RigidBody, 2) ShearCell
Ini         =   (T=:lineara,)
# ------------------------------------------------------------------- #
# Plot Einstellungen ================================================ #
Pl  =   (
    qinc        =   5,
    qsc         =   1.0 # 2e2*(100*(60*60*24*365.15)),
)
# Animationssettings ================================================ #
k           =   scatter()
path        =   string("./Results/")
anim        =   Plots.Animation(path, String[] )
save_fig    =   0
# ------------------------------------------------------------------- #
# Modellgeometrie Konstanten ======================================== #
M   =   (
    xmin    =   [0.0],              #   [ m ] 
    xmax    =   [8700e3],           #   [ m ]
    ymin    =   [-2900e3],          #   [ m ]
    ymax    =   [0.0],              #   [ m ]
)
# ------------------------------------------------------------------- #
# Referenzparameter ================================================= #
P   =   (
    g   =   9.81,                   #   Schwerebeschleunigung [m/s^2]
    ρ₀  =   3300.0,                 #   Hintergunddichte [kg/m^3]
    k   =   4.125,                  #   Thermische Leitfaehigkeit [ W/m/K ]
    cp  =   1250.0,                 #   Heat capacity [ J/kg/K ]
    α   =   2.0e-5,                 #   Thermischer Expnasionskoef. [ K^-1 ]
    Q₀  =   [0.0],                  #   Waermeproduktionsrate pro Volumen [W/m^3]
    η₀  =   [3.947725485e23],       #   Viskositaet [ Pa*s ] [1.778087025e21]
)
P1  =   (
    κ       =   P.k/P.ρ₀/P.cp,      # 	Thermische Diffusivitaet [ m^2/s ]
    ΔT      =   2500.0,             #   Temperaturdifferenz
    # Falls Ra < 0 gesetzt ist, dann wird Ra aus den obigen Parametern
    # berechnet. Falls Ra gegeben ist, dann wird die Referenzviskositaet so
    # angepasst, dass die Skalierungsparameter die gegebene Rayleigh-Zahl
    # ergeben.
    Ra      =   [1.0e4],       #   Rayleigh number
)
P2  =   (
    Ttop    =   [273.15],           #   Temperatur an der Oberfläche [ K ]
    Tbot    =   [273.15 + P1.ΔT],   #   Temperatur an der Unterseite [ K ]
)
P   =   merge(P,P1,P2)
filename    =   string("11_ThermalConvection_",P.Ra[1],
                        "_",Ini.T,"_",FD.Method.Adv,"_",FD.Method.Diff,
                        "_",FD.Method.Mom)
# ------------------------------------------------------------------- #
# Scaling =========================================================== # 
# Scaling constants ---
S   =   (
    hsc =   (M.ymax[1]-M.ymin[1]),                          #   Length scale [ m ]
    tsc =   (M.ymax[1]-M.ymin[1])^2 / P.κ,                  #   Time scale [ s ]
    vsc =   P.κ / (M.ymax[1]-M.ymin[1]),                    #   Velocity scale [ m/s ]
    τsc =   (P.η₀[1] * P.κ)/(M.ymax[1]-M.ymin[1]),          #   Stress scale [ Pa ]
    Tsc =   P.ΔT,                                           #   Temperature scale [ K ]
    Qsc =   (P.ΔT*P.κ*P.ρ₀*P.cp)/(M.ymax[1]-M.ymin[1])^2,   #   Heat source scale [ w/m³ ]
)
# ------------------------------------------------------------------- #
# Grid ============================================================== #
NC  =   (
    x   =   150,    #   horizontal grid resolution
    y   =   50,    #   vertical grid resolution
)
NV      =   (
    x   =   NC.x + 1,
    y   =   NC.y + 1,
)
Δ       =   (
    x   =   [(M.xmax[1] - M.xmin[1])/NC.x],
    y   =   [(M.ymax[1] - M.ymin[1])/NC.y],
)
# ------------------------------------------------------------------- #
# Allocation ======================================================== #
D       =   (
    Q       =   zeros(Float64,(NC...)),
    T       =   zeros(Float64,(NC...)),
    T0      =   zeros(Float64,(NC...)),
    T_ex    =   zeros(Float64,(NC.x+2,NC.y+2)),
    T_exo   =   zeros(Float64,(NC.x+2,NC.y+2)),
    ρ       =   ones(Float64,(NC...)),
    cp      =   ones(Float64,(NC...)),
    vx      =   zeros(Float64,(NV.x,NV.y+1)),
    vy      =   zeros(Float64,(NV.x+1,NV.y)),    
    Pt      =   zeros(Float64,(NC...)),
    vxc     =   zeros(Float64,(NC...)),
    vyc     =   zeros(Float64,(NC...)),
    vc      =   zeros(Float64,(NC...)),
    wt      =   zeros(Float64,(NC...)),
    wtv     =   zeros(Float64,(NV...)),
    ΔTtop   =   zeros(Float64,NC.x),
    ΔTbot   =   zeros(Float64,NC.x),
    Tmax    =   [0.0],
    Tmin    =   [0.0],
    Tmean   =   [0.0],
)
# Heat production rate ------
@. D.Q      =   P.Q₀
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
# Residuals ---
Fm     =    (
    x       =   zeros(Float64,NV.x, NC.y), 
    y       =   zeros(Float64,NC.x, NV.y)
)
FPt         =   zeros(Float64,NC...)
# ------------------------------------------------------------------- #
# Time ============================================================== #
T   =   (
    year    =   365.25*3600*24,     #   Seconds per year
    tmax    =   [1000000.0],           #   [ Ma ]
    Δfacc   =   0.9,                #   Courant time factor
    Δfacd   =   0.9,                #   Diffusion time factor
    Δ       =   [0.0],
    Δc      =   [0.0],              #   Courant time step
    Δd      =   [0.0],              #   Diffusion time stability criterion
    itmax   =   8000,              #   Maximum iterations; 30000
)
T.tmax[1]   =   T.tmax[1]*1e6*T.year    #   [ s ]
T.Δc[1]     =   T.Δfacc * minimum((Δ.x[1],Δ.y[1])) / 
                    (sqrt(maximum(abs.(D.vx))^2 + maximum(abs.(D.vy))^2))
T.Δd[1]     =   T.Δfacd * (1.0 / (2.0 * P.κ *(1.0/Δ.x[1]^2 + 1/Δ.y[1]^2)))

T.Δ[1]      =   minimum([T.Δd[1],T.Δc[1]])

Time        =   zeros(T.itmax)
Nus         =   zeros(T.itmax)
meanV       =   zeros(T.itmax)
meanT       =   zeros(T.itmax,NC.y+2)
find        =   0
# ------------------------------------------------------------------- #
# Scaling laws ====================================================== #
# Geometry ---
M.xmin[1]   /=  S.hsc
M.xmax[1]   /=  S.hsc
M.ymin[1]   /=  S.hsc
M.ymax[1]   /=  S.hsc
Δ.x[1]      /=  S.hsc
Δ.y[1]      /=  S.hsc         
# Time --- 
T.tmax[1]   /=  S.tsc
T.Δc[1]     /=  S.tsc
T.Δd[1]     /=  S.tsc
T.Δ[1]      /=  S.tsc
# Temperature etc. ---
P.Ttop[1]   =   (P.Ttop[1] - 273.15)/S.Tsc
P.Tbot[1]   =   (P.Tbot[1] - 273.15)/S.Tsc
@. D.Q      /=  S.Qsc
# # Viscosity ---
# @. D.η      /=  P.η₀[1]
# ------------------------------------------------------------------- #
# Coordinates ======================================================= #
x       =   (
    c   =   LinRange(M.xmin[1]+Δ.x[1]/2,M.xmax[1]-Δ.x[1]/2,NC.x),
    ce  =   LinRange(M.xmin[1] - Δ.x[1]/2.0, M.xmax[1] + Δ.x[1]/2.0, NC.x+2),
    v   =   LinRange(M.xmin[1],M.xmax[1],NV.x),
)
y       =   (
    c   =   LinRange(M.ymin[1]+Δ.y[1]/2,M.ymax[1]-Δ.y[1]/2,NC.y),
    ce  =   LinRange(M.ymin[1] - Δ.x[1]/2.0, M.ymax[1] + Δ.x[1]/2.0, NC.y+2),
    v   =   LinRange(M.ymin[1],M.ymax[1],NV.y),
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
# ------------------------------------------------------------------- #
# Anfangsbedingungen ================================================ #
# Temperatur ------
# Ttop    =   273.15
# Tbot    =   Ttop + P.ΔT
if FD.Method.Adv==:tracers 
    # Tracer Initialization ---
    # Need to implement incremental marker update first! 
    nmx,nmy     =   3,3
    noise       =   0
    nmark       =   nmx*nmy*NC.x*NC.y
    Aparam      =   :thermal
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
    Ma      =   IniTracer2D(Aparam,nmx,nmy,Δ,M,NC,noise,Ini.p,phase)
    # RK4 weights ---
    rkw     =   1.0/6.0*[1.0 2.0 2.0 1.0]   # for averaging
    rkv     =   1.0/2.0*[1.0 1.0 2.0 2.0]   # for time stepping
    # Count marker per cell ---
    CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,NV,1)
    # Interpolate from markers to cell ---
    Markers2Cells(Ma,nmark,MPC.PG_th,D.ρ,MPC.wt_th,D.wt,x,y,Δ,Aparam,ρ)
else
    IniTemperature!(Ini.T,M,NC,Δ,D,x,y;Tb=P.Tbot[1],Ta=P.Ttop[1])
    if FD.Method.Adv==:slf
        D.T_exo    .=  D.T_ex
    end
end
# ------------------------------------------------------------------- #
# Boundary Conditions =============================================== #
# Temperature ------
TBC     = (
    type    = (W=:Neumann, E=:Neumann,N=:Dirichlet,S=:Dirichlet),
    val     = (W=zeros(NC.y),E=zeros(NC.y),
                    N=P.Ttop[1].*ones(NC.x),S=P.Tbot[1].*ones(NC.x)))
# Velocity ------
VBC     =   (
    type    =   (E=:freeslip,W=:freeslip,S=:freeslip,N=:freeslip),
    val     =   (E=zeros(NV.y),W=zeros(NV.y),S=zeros(NV.x),N=zeros(NV.x)),
)
# ------------------------------------------------------------------- # 
# Rayleigh Zahl Bedingungen ========================================= #
if P.Ra[1] < 0
    # Falls die Rayleigh Zahl nicht explizit angegeben wird, dann 
    # wird sie hier berechnet
    P.Ra[1]     =   P.ρ₀*P.g*P.α*P.ΔT*S.hsc^3/P.η₀[1]/P.κ
else
    # Falls die Rayleigh Zahl explizit angegeben ist, dann wird hier 
    # die Referenzviskositaet η₀ angepasst. 
    P.η₀[1]     =   P.ρ₀*P.g*P.α*P.ΔT*S.hsc^3/P.Ra[1]/P.κ
end
# =================================================================== #
# Linear Equations ================================================== #
# Momentum Equation ======
off    = [  NV.x*NC.y,                          # vx
            NV.x*NC.y + NC.x*NV.y,              # vy
            NV.x*NC.y + NC.x*NV.y + NC.x*NC.y]  # Pt

Num    =    (
    Vx  =   reshape(1:NV.x*NC.y, NV.x, NC.y), 
    Vy  =   reshape(off[1]+1:off[1]+NC.x*NV.y, NC.x, NV.y), 
    Pt  =   reshape(off[2]+1:off[2]+NC.x*NC.y,NC...),
    T   =   reshape(1:NC.x*NC.y, NC.x, NC.y),
)
ndof    =   maximum(Num.T)        
# Temperature Equation ======
if FD.Method.Diff==:implicit || FD.Method.Diff==:CNA
    if FD.Method.Diff==:CNA
        K1      =   ExtendableSparseMatrix(ndof,ndof)
        K2      =   ExtendableSparseMatrix(ndof,ndof)
    else
        K       =   ExtendableSparseMatrix(ndof,ndof)
    end
    rhs         =   zeros(ndof)
elseif FD.Method.Diff==:dc
    niter       =   10
    ϵ           =   1e-10
    k           =   (x=ones(NC.x+1,NC.y), y=ones(NC.x,NC.y+1))
    K           =   ExtendableSparseMatrix(ndof,ndof)
    R           =   zeros(Float64,NC...)
    ∂T          =   (∂x=zeros(NC.x+1, NC.y), ∂y=zeros(NC.x, NC.y+1))
    q           =   (x=zeros(NC.x+1, NC.y), y=zeros(NC.x, NC.y+1))
end
# ------------------------------------------------------------------- #
# Time Loop ========================================================= #
for it = 1:T.itmax
    χ       =   zeros(maximum(Num.Pt))      #   Unknown Vector ME
    rhsM    =   zeros(maximum(Num.Pt))      #   Right-hand Side ME
    if it>1
        Time[it]  =   Time[it-1] + T.Δ[1]
    end
    @printf("Time step: #%04d, Time [non-dim]: %04e\n ",it,
                    Time[it])
    # Momentum Equation(ME) =======
    D.vx    .=  0.0
    D.vy    .=  0.0 
    D.Pt    .=  0.0
    if FD.Method.Mom==:direct
        # Update K ---
        KM      =   Assemblyc(NC, NV, Δ, 1.0, VBC, Num)
        # Update RHS ---
        # rhs term defined by the Boussinesq approximation
        rhsM    =   updaterhsc( NC, NV, Δ, 1.0, -P.Ra[1]*D.T, -1.0, VBC, Num )
        # Solve System of Equations ---
        χ       =   KM \ rhsM
        # Update Unknown Variables ---
        D.vx[:,2:end-1]     .=  χ[Num.Vx]
        D.vy[2:end-1,:]     .=  χ[Num.Vy]
        D.Pt                .=  χ[Num.Pt]
    elseif FD.Method.Mom==:dc
        # Initial Residual -------------------------------------------------- #
        @. D.ρ  =   -P.Ra[1]*D.T
        Residuals2Dc!(D,VBC,ε,τ,divV,Δ,1.0,1.0,Fm,FPt)
        rhsM[Num.Vx]    =   Fm.x[:]
        rhsM[Num.Vy]    =   Fm.y[:]
        rhsM[Num.Pt]    =   FPt[:]
        # ------------------------------------------------------------------- #
        # Assemble Coefficients ============================================= #
        K       =   Assemblyc(NC, NV, Δ, 1.0, VBC, Num)
        # ------------------------------------------------------------------- #
        # Solution of the linear system ===================================== #
        χ      =   - K \ rhsM
        # ------------------------------------------------------------------- #
        # Update Unknown Variables ========================================== #
        D.vx[:,2:end-1]     .+=  χ[Num.Vx]
        D.vy[2:end-1,:]     .+=  χ[Num.Vy]
        D.Pt                .+=  χ[Num.Pt]
        @. D.ρ  =   ones(NC...)
    end
    # ======
    # Get the velocity on the centroids ------
    for i = 1:NC.x
        for j = 1:NC.y
            D.vxc[i,j]  = (D.vx[i,j+1] + D.vx[i+1,j+1])/2
            D.vyc[i,j]  = (D.vy[i+1,j] + D.vy[i+1,j+1])/2
        end
    end
    @. D.vc        = sqrt(D.vxc^2 + D.vyc^2)
    # ---
    @show(maximum(D.vc))
    @show(minimum(D.Pt))
    @show(maximum(D.Pt))
    # Calculate time stepping ======================================= #
    T.Δc[1]     =   T.Δfacc * minimum((Δ.x[1],Δ.y[1])) / 
            (sqrt(maximum(abs.(D.vx))^2 + maximum(abs.(D.vy))^2))
    T.Δd[1]     =   T.Δfacd * (1.0 / (2.0 *(1.0/Δ.x[1]^2 + 1/Δ.y[1]^2)))
    T.Δ[1]      =   minimum([T.Δd[1],T.Δc[1]])
    if Time[it] > T.tmax[1] 
        T.Δ[1]      =   T.tmax[1] - Time[it-1]
        Time[it]    =   Time[it-1] + T.Δ[1]
        it          =   T.itmax
    end
    # Plot ========================================================== #
    if mod(it,10) == 0 || it == T.itmax || it == 1
        p = heatmap(x.c,y.c,D.T',
                xlabel="x",ylabel="y",colorbar=true,
                title="Temperature",color=cgrad(:lajolla),
                aspect_ratio=:equal,xlims=(M.xmin[1], M.xmax[1]),
                ylims=(M.ymin[1], M.ymax[1]),
                layout=(2,1),subplot=1)
        heatmap!(p,x.c,y.c,D.vc',color=:imola,
            xlabel="x[km]",ylabel="y[km]",colorbar=true,
            title="Velocity",aspect_ratio=:equal,
            xlims=(M.xmin[1], M.xmax[1]), 
            ylims=(M.ymin[1], M.ymax[1]),
            layout=(2,1),subplot=2)
        quiver!(p,x.c2d[1:Pl.qinc:end,1:Pl.qinc:end],
            y.c2d[1:Pl.qinc:end,1:Pl.qinc:end],
            quiver=(D.vxc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc,
                    D.vyc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc),
            la=0.5,color="black",
            layout=(2,1),subplot=2)
        if save_fig == 1
            Plots.frame(anim)
        elseif save_fig == 0
            display(p)
        end
    end
    # --------------------------------------------------------------- #
    # Advection ===================================================== #
    if FD.Method.Adv==:upwind
        upwindc2D!(D.T,D.T_ex,D.vxc,D.vyc,NC,T.Δ[1],Δ.x[1],Δ.y[1])            
    elseif FD.Method.Adv==:slf
        slfc2D!(D.T,D.T_ex,D.T_exo,D.vxc,D.vyc,NC,T.Δ[1],Δ.x[1],Δ.y[1])
    elseif FD.Method.Adv==:semilag
        semilagc2D!(D.T,D.T_ex,D.vxc,D.vyc,[],[],x,y,T.Δ[1])
    elseif FD.Method.Adv==:tracers
        # Advect tracers ---
        @printf("Running on %d thread(s)\n", nthreads())  
        AdvectTracer2D(Ma,nmark,D,x,y,T.Δ[1],Δ,NC,rkw,rkv,1)
        CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,NV,it)
        # Interpolate phase from tracers to grid ---
        Markers2Cells(Ma,nmark,MPC.PG_th,D.ρ,MPC.wt_th,D.wt,x,y,Δ,Aparam,ρ)
    end
    # --------------------------------------------------------------- #
    # Diffusion ===================================================== #
    if FD.Method.Diff==:explicit
        ForwardEuler2Dc!(D, 1.0, Δ.x[1], Δ.y[1], T.Δ[1], D.ρ, 1.0, NC, TBC)
    elseif FD.Method.Diff==:implicit
        BackwardEuler2Dc!(D, 1.0, Δ.x[1], Δ.y[1], T.Δ[1], D.ρ, 1.0, NC, TBC, rhs, K, Num)
    elseif FD.Method.Diff==:CNA
        CNA2Dc!(D, 1.0, Δ.x[1], Δ.y[1], T.Δ[1], D.ρ, 1.0, NC, TBC, rhs, K1, K2, Num)
    elseif FD.Method.Diff==:ADI
        ADI2Dc!(D, 1.0, Δ.x[1], Δ.y[1], T.Δ[1], D.ρ, 1.0, NC, TBC)
    elseif FD.Method.Diff==:dc
        D.T0    .=  D.T
        for iter = 1:niter
            # Evaluate residual
            ComputeResiduals2D!(R, D.T, D.T_ex, D.T0, ∂T, q, D.ρ, D.cp, k, TBC, Δ, T.Δ[1])
            # @printf("||R|| = %1.4e\n", norm(R)/length(R))
            norm(R)/length(R) < ϵ ? break : nothing
            # Assemble linear system
            K  = AssembleMatrix2D(D.ρ, D.cp, k, TBC, Num, NC, Δ, T.Δ[1])
            # Solve for temperature correction: Cholesky factorisation
            Kc = cholesky(K.cscmatrix)
            # Solve for temperature correction: Back substitutions
            δT = -(Kc\R[:])
            # Update temperature
            @. D.T += δT[Num.T]
        end        
    end
    # --------------------------------------------------------------- #
    @printf("\n")
    # Heat flow at the surface ====================================== #
    @. D.ΔTbot  =   
         (((D.T_ex[2:end-1,2]+D.T_ex[2:end-1,3])/2.0) - 
        ((D.T_ex[2:end-1,2]+D.T_ex[2:end-1,1])/2.0)) / Δ.y[1]
    @. D.ΔTtop  =   
        (((D.T_ex[2:end-1,end-2]+D.T_ex[2:end-1,end-1]) / 2.0) - 
        ((D.T_ex[2:end-1,end-1]+D.T_ex[2:end-1,end]) / 2.0)) / Δ.y[1]
    Nus[it]     =   mean(D.ΔTtop)
    meanT[it,:] =   mean(D.T_ex,dims=1)
    meanV[it]   =   mean(D.vc)
    # --------------------------------------------------------------- #
    # Check break =================================================== #
    # If the maximum time is reached or if the models reaches steady, 
    # state the time loop is stoped! 
    if Time[it] > 0.001
        epsC    =   1e-16; 
        ind     =   findfirst(Time .> 
                        (Time[it] - 0.001))
        epsV    =   std(meanV[ind:it])
        # @show epsV1, epsV
        if save_fig == 1
            @show it,log10((epsV))
            plot!(k,(it,log10((epsV))),
                xlabel="it",ylabel="log₁₀(εᵥ)",label="",
                markershape=:circle,markercolor=:black)
        end
        find    =   it
        @printf("ε_V = %g, ε_C = %g \n",epsV,epsC)
        if Time[it] >= T.tmax[1]
            @printf("Maximum time reached!\n")
            find    =   it
            break
        elseif (epsV <= epsC)
            @printf("Convection reaches steady state!\n")
            find    =   it
            break
        end
    end
    # --------------------------------------------------------------- #
    @printf("\n")
end
if save_fig == 1
    # Write the frames to a GIF file
    Plots.gif(anim, string( path, filename, ".gif" ), fps = 15)
    foreach(rm, filter(startswith(string(path,"00")), readdir(path,join=true)))
end
# Save final figure ===================================================== #
p2 = heatmap(x.c,y.c,D.T',
            xlabel="x[km]",ylabel="y[km]",colorbar=true,
            title="Temperature",color=cgrad(:lajolla),
            aspect_ratio=:equal,xlims=(M.xmin[1], M.xmax[1]),
            ylims=(M.ymin[1], M.ymax[1]),
            layout=(2,1),subplot=1)
    heatmap!(p2,x.c,y.c,D.vc',color=:imola,
            xlabel="x[km]",ylabel="y[km]",colorbar=true,
            title="Velocity",aspect_ratio=:equal,
            xlims=(M.xmin[1], M.xmax[1]), 
            ylims=(M.ymin[1], M.ymax[1]),
            layout=(2,1),subplot=2)
    quiver!(p2,x.c2d[1:Pl.qinc:end,1:Pl.qinc:end],
            y.c2d[1:Pl.qinc:end,1:Pl.qinc:end],
            quiver=(D.vxc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc,
                    D.vyc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc),
            la=0.5,color="black",
            layout=(2,1),subplot=2)
if save_fig == 1
    savefig(k,string("./Results/11_ThermalConvection_iterations_",P.Ra[1],
            "_",Ini.T,"_",FD.Method.Adv,"_",FD.Method.Diff,"_",FD.Method.Mom,".png"))
    savefig(p2,string("./Results/11_ThermalConvection_Final_Stage_",P.Ra[1],"_it_",find,"_",
                        Ini.T,"_",FD.Method.Adv,"_",FD.Method.Diff,"_",FD.Method.Mom,".png"))
elseif save_fig == 0
    display(p2)
end
# ----------------------------------------------------------------------- #
# Plot time serieses ==================================================== #
q2  =   plot(Time[1:find],Nus[1:find],
            xlabel="Time [ non-dim ]", ylabel="Nus",label="",
            layout=(2,1),suplot=1)
plot!(q2,Time[1:find],meanV[1:find],
            xlabel="Time [ non-dim ]", ylabel="V_{RMS}",label="",
            layout=(2,1),subplot=2)
if save_fig == 1
    savefig(q2,string("./Results/11_ThermalConvectionTimeSeries",P.Ra[1],"_",
                        Ini.T,"_",FD.Method.Adv,"_",FD.Method.Diff,"_",FD.Method.Mom,".png"))
elseif save_fig == 0
    display(q2)
end
# ======================================================================= #
end

ThermalConvection_scaled()
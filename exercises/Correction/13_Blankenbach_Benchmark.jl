using Plots, ExtendableSparse
using GeoModBox
using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
using GeoModBox.AdvectionEquation.TwoD, GeoModBox.HeatEquation.TwoD
using GeoModBox.Scaling
using Statistics, LinearAlgebra
using Printf

function BlankenbachBenchmark()
# Define Initial Condition ============================================== #
# Temperature - 
#   1) circle, 2) gaussian, 3) block, 4) linear, 5) lineara
# !!! Gaussian is not working!!! 
Ini         =   (T=:lineara,)
# ----------------------------------------------------------------------- #
# Plot Settings ========================================================= #
Pl  =   (
    qinc        =   5,
    qsc         =   1.0e-3
)
# Animation Settings ==================================================== #
k           =   scatter()
path        =   string("./exercises/Correction/Results/")
anim        =   Plots.Animation(path, String[] )
save_fig    =   1
# ----------------------------------------------------------------------- #
# Geometry ============================================================== #
M   =   Geometry(
    xmin    =   0.0,                #   [ m ] 
    xmax    =   1000e3,             #   [ m ]
    ymin    =   -1000e3,            #   [ m ]
    ymax    =   0.0,                #   [ m ]
)
# ----------------------------------------------------------------------- #
# Referenzparameter ===================================================== #
P   =   Physics(
    g       =   10.0,               #   Graviational Acceleration [m/s^2]
    ρ₀      =   4000.0,             #   Reference Density [kg/m^3]
    k       =   5.0,                #   Thermal Conductivity [ W/m/K ]
    cp      =   1250.0,             #   Heat capacity [ J/kg/K ]
    α       =   2.5e-5,             #   THermal Expansion [ K^-1 ]
    Q₀      =   0.0,                #   Volumetric Heat Production rate [W/m^3]
    η₀      =   1e21,               #   Reference Viscosity [ Pa*s ]
    ΔT      =   1000.0,             #   Temperature Difference
    # If Ra < 0, Ra will be calculated from the reference parameters.
    # If Ra is defined, the reference viscosity will be adjusted such that
    # the scaling parameters result in the given Ra.
    Ra      =   -9999,              #   Rayleigh number
    Ttop    =   273.15,             #   Temperature at the surface [ K ]
)
# ----------------------------------------------------------------------- #
# Define Scaling Constants ============================================== # 
S   =   ScalingConstants!(M,P)
# ----------------------------------------------------------------------- #
# Numerical Grid ======================================================== #
NC  =   (
    x   =   50,
    y   =   50,
)
NV      =   (
    x   =   NC.x + 1,
    y   =   NC.y + 1,
)
Δ       =   GridSpacing(
    x   =   (M.xmax - M.xmin)/NC.x,
    y   =   (M.ymax - M.ymin)/NC.y,
)
# ----------------------------------------------------------------------- #
# Data Arrays =========================================================== #
D       =   DataFields(
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
    Tmax    =   0.0,
    Tmin    =   0.0,
    Tmean   =   0.0,
)
# Heat Production Rate ------
@. D.Q      =   P.Q₀
# ----------------------------------------------------------------------- #
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
# ----------------------------------------------------------------------- #
# Time parameters ======================================================= #
T   =   TimeParameter(
    tmax    =   1000000.0,          #   [ Ma ]
    Δfacc   =   0.9,                #   Courant time factor
    Δfacd   =   0.9,                #   Diffusion time factor
    itmax   =   8000,               #   Maximum iterations; 30000
)
T.tmax      =   T.tmax*1e6*T.year   #   [ s ]
T.Δc        =   T.Δfacc * minimum((Δ.x,Δ.y)) / 
                    (sqrt(maximum(abs.(D.vx))^2 + maximum(abs.(D.vy))^2))
T.Δd        =   T.Δfacd * (1.0 / (2.0 * P.κ *(1.0/Δ.x^2 + 1/Δ.y^2)))

T.Δ         =   minimum([T.Δd,T.Δc])

Time        =   zeros(T.itmax)
Nus         =   zeros(T.itmax)
meanV       =   zeros(T.itmax)
meanT       =   zeros(T.itmax,NC.y+2)
find        =   0
# ----------------------------------------------------------------------- #
# Scaling laws ========================================================== #
ScaleParameters!(S,M,Δ,T,P,D)
# ----------------------------------------------------------------------- #
# Coordinates =========================================================== #
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
# Initial Condition ===================================================== #
# Temperature ------
IniTemperature!(Ini.T,M,NC,D,x,y;Tb=P.Tbot,Ta=P.Ttop)
# ----------------------------------------------------------------------- #
# Boundary Conditions =================================================== #
# Temperature ------
TBC     = (
    type    = (W=:Neumann, E=:Neumann,N=:Dirichlet,S=:Dirichlet),
    val     = (W=zeros(NC.y),E=zeros(NC.y),
                    N=P.Ttop.*ones(NC.x),S=P.Tbot.*ones(NC.x)))
# Velocity ------
VBC     =   (
    type    =   (E=:freeslip,W=:freeslip,S=:freeslip,N=:freeslip),
    val     =   (E=zeros(NV.y),W=zeros(NV.y),S=zeros(NV.x),N=zeros(NV.x)),
)
# ----------------------------------------------------------------------- #
# Rayleigh Number ======================================================= #
if P.Ra < 0
    P.Ra     =   P.ρ₀*P.g*P.α*P.ΔT*S.hsc^3/P.η₀/P.κ
else
    P.η₀     =   P.ρ₀*P.g*P.α*P.ΔT*S.hsc^3/P.Ra/P.κ
end
filename    =   string("13_Blankenbach_",@sprintf("%.2e",P.Ra),
                        "_",NC.x,"_",NC.y,
                        "_",Ini.T)
# ----------------------------------------------------------------------- #
# Linear System of Equations ============================================ #
# Momentum Conservation Equation (MCE) ------
niter  =   50
ϵ      =   1e-10
off    = [  NV.x*NC.y,                          # vx
            NV.x*NC.y + NC.x*NV.y,              # vy
            NV.x*NC.y + NC.x*NV.y + NC.x*NC.y ] # Pt
Num    =    (
    Vx  =   reshape(1:NV.x*NC.y, NV.x, NC.y), 
    Vy  =   reshape(off[1]+1:off[1]+NC.x*NV.y, NC.x, NV.y), 
    Pt  =   reshape(off[2]+1:off[2]+NC.x*NC.y,NC...),
    T   =   reshape(1:NC.x*NC.y, NC.x, NC.y),
)
ndof    =   maximum(Num.T)        
# Energy Conservation Equation (ECE) ------
K1      =   ExtendableSparseMatrix(ndof,ndof)
K2      =   ExtendableSparseMatrix(ndof,ndof)
rhs     =   zeros(ndof)
# ----------------------------------------------------------------------- #
# Time Loop ============================================================= #
for it = 1:T.itmax
    χ       =   zeros(maximum(Num.Pt))      #   Unknown Vector MCE
    rhsM    =   zeros(maximum(Num.Pt))      #   Right Hand Side MCE
    if it>1
        Time[it]  =   Time[it-1] + T.Δ
    end
    @printf("Time step: #%04d, Time [non-dim]: %04e\n ",it,
                    Time[it])
    # MCE ------
    D.vx    .=  0.0
    D.vy    .=  0.0 
    D.Pt    .=  0.0
    # Residual Calculation ------
    @. D.ρ  =   -P.Ra*D.T
    for iter = 1:niter
        Residuals2Dc!(D,VBC,ε,τ,divV,Δ,1.0,1.0,Fm,FPt)
        rhsM[Num.Vx]    =   Fm.x[:]
        rhsM[Num.Vy]    =   Fm.y[:]
        rhsM[Num.Pt]    =   FPt[:]
        @printf("||R_M|| = %1.4e\n", norm(rhsM)/length(rhsM))
            norm(rhsM)/length(rhsM) < ϵ ? break : nothing
        # Update K ------
        K       =   Assemblyc(NC, NV, Δ, 1.0, VBC, Num)
        # Solving Linear System of Equations ------
        χ      =   - K \ rhsM
        # Update unknown variables ------
        D.vx[:,2:end-1]     .+=  χ[Num.Vx]
        D.vy[2:end-1,:]     .+=  χ[Num.Vy]
        D.Pt                .+=  χ[Num.Pt]
    end
    D.ρ  =   ones(NC...)
    # ======
    # Calculate Velocity on the Centroids ------
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
    # Calculate time stepping =========================================== #
    T.Δc        =   T.Δfacc * minimum((Δ.x,Δ.y)) / 
            (sqrt(maximum(abs.(D.vx))^2 + maximum(abs.(D.vy))^2))
    T.Δd        =   T.Δfacd * (1.0 / (2.0 *(1.0/Δ.x^2 + 1/Δ.y^2)))
    T.Δ         =   minimum([T.Δd,T.Δc])
    if Time[it] > T.tmax
        T.Δ         =   T.tmax - Time[it-1]
        Time[it]    =   Time[it-1] + T.Δ
        it          =   T.itmax
    end
    # Plot ============================================================== #
    if mod(it,10) == 0 || it == T.itmax || it == 1
        p = heatmap(x.c,y.c,D.T',
                xlabel="x",ylabel="y",colorbar=true,
                title="Temperature",color=cgrad(:lajolla),
                aspect_ratio=:equal,xlims=(M.xmin, M.xmax),
                ylims=(M.ymin, M.ymax),
                layout=(2,1),subplot=1)
        heatmap!(p,x.c,y.c,D.vc',color=:imola,
            xlabel="x[km]",ylabel="y[km]",colorbar=true,
            title="Velocity",aspect_ratio=:equal,
            xlims=(M.xmin, M.xmax), 
            ylims=(M.ymin, M.ymax),
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
    # ------------------------------------------------------------------- #
    # Advection ========================================================= #
    semilagc2D!(D.T,D.T_ex,D.vxc,D.vyc,[],[],x,y,T.Δ)
    # ------------------------------------------------------------------- #
    # Diffusion ========================================================= #
    CNA2Dc!(D, 1.0, Δ.x, Δ.y, T.Δ, D.ρ, 1.0, NC, TBC, rhs, K1, K2, Num)
    # ------------------------------------------------------------------- #
    # Surface Heat Flow ================================================= #
    @. D.ΔTbot  =   
         (((D.T_ex[2:end-1,2]+D.T_ex[2:end-1,3])/2.0) - 
        ((D.T_ex[2:end-1,2]+D.T_ex[2:end-1,1])/2.0)) / Δ.y
    @. D.ΔTtop  =   
        (((D.T_ex[2:end-1,end-2]+D.T_ex[2:end-1,end-1]) / 2.0) - 
        ((D.T_ex[2:end-1,end-1]+D.T_ex[2:end-1,end]) / 2.0)) / Δ.y
    Nus[it]     =   mean(D.ΔTtop)
    meanT[it,:] =   mean(D.T_ex,dims=1)
    meanV[it]   =   mean(D.vc)
    # ------------------------------------------------------------------- #
    # Check break ======================================================= #
    # If the maximum time is reached or if the models reaches steady 
    # state the time loop is stoped! 
    if Time[it] > 0.0038
        epsC    =   1e-3; 
        ind     =   findfirst(Time .> 
                        (Time[it] - 0.0038))
        epsV    =   std(meanV[ind:it])
        if save_fig == 1
            plot!(k,(it,log10((epsV))),
                xlabel="it",ylabel="log₁₀(εᵥ)",label="",
                markershape=:circle,markercolor=:black)
        end
        find    =   it
        @printf("ε_V = %g, ε_C = %g \n",epsV,epsC)
        if Time[it] >= T.tmax
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
            aspect_ratio=:equal,xlims=(M.xmin, M.xmax),
            ylims=(M.ymin, M.ymax),
            layout=(2,1),subplot=1)
    heatmap!(p2,x.c,y.c,D.vc',color=:imola,
            xlabel="x[km]",ylabel="y[km]",colorbar=true,
            title="Velocity",aspect_ratio=:equal,
            xlims=(M.xmin, M.xmax), 
            ylims=(M.ymin, M.ymax),
            layout=(2,1),subplot=2)
    quiver!(p2,x.c2d[1:Pl.qinc:end,1:Pl.qinc:end],
            y.c2d[1:Pl.qinc:end,1:Pl.qinc:end],
            quiver=(D.vxc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc,
                    D.vyc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc),
            la=0.5,color="black",
            layout=(2,1),subplot=2)
if save_fig == 1
    savefig(k,string("./exercises/Correction/Results/13_BlankenbachBenchmark_Scaled_iterations_",@sprintf("%.2e",P.Ra),
            "_",NC.x,"_",NC.y,
            "_",Ini.T,".png"))
    savefig(p2,string("./exercises/Correction/Results/13_BlankenbachBenchmark_Scaled_Final_Stage_",@sprintf("%.2e",P.Ra),
            "_",NC.x,"_",NC.y,"_it_",find,"_",
            Ini.T,".png"))
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
    savefig(q2,string("./exercises/Correction/Results/13_BlankenbachBenchmarkTimeSeries_Scaled_",@sprintf("%.2e",P.Ra),
                        "_",NC.x,"_",NC.y,"_",Ini.T,".png"))
elseif save_fig == 0
    display(q2)
end
# ======================================================================= #
end

BlankenbachBenchmark()
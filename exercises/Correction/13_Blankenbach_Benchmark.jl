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
    qinc        =   10,
    qsc         =   5.0e-5 # 1.0e-3, 4.0e-4, 5.0e-5
)
# Animation Settings ==================================================== #
k           =   scatter()
path        =   string("./exercises/Correction/Results/")
anim        =   Plots.Animation(path, String[] )
save_fig    =   0
# ----------------------------------------------------------------------- #
# Benchmark Values ====================================================== # 
# Taken from Gerya (2019), Introduction to numerical geodynamic modelling
B       =   (
    Nr      =   3,
    # Nusselt Number at the top
    Nu      =   [4.8844,10.534,21.972,10.066,6.9229],   
    # Root Mean Square Velocity-scaled
    Vrms    =   [42.865,193.21,833.99,480.43,171.76],   
    # Non-dimensional temperature gradients in the model corner
    q1      =   [8.0593,19.079,45.964,17.531,18.484],
    q2      =   [0.5888,0.7228,0.8772,1.0085,0.1774],
    q3      =   [8.0593,19.079,45.964,26.809,14.168], 
    q4      =   [0.5888,0.7228,0.8772,0.4974,0.6177],
    # Local minimum along the central vertical temperature profile 
    Tmin    =   [0.4222,0.4284,0.4322,0.7405,0.3970],
    ymin    =   [0.2249,0.1118,0.0577,0.0623,0.1906],
    # Local maximum algon the central vertical temperature profile
    Tmax    =   [0.5778,0.5716,0.5678,0.8323,0.5758],
    ymax    =   [0.7751,0.8882,0.9423,0.8243,0.7837],
)
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
    Ttop    =   273,                #   Temperature at the surface [ K ]
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
    T_ex0   =   zeros(Float64,(NC.x+2,NC.y+2)),
    Told_ex =   zeros(Float64,(NC.x+2,NC.y+2)),
    ρ       =   ones(Float64,(NC...)),
    vx      =   zeros(Float64,(NV.x,NV.y+1)),
    vy      =   zeros(Float64,(NV.x+1,NV.y)),    
    Pt      =   zeros(Float64,(NC...)),
    vxc     =   zeros(Float64,(NC...)),
    vyc     =   zeros(Float64,(NC...)),
    vxco    =   zeros(Float64,(NC...)),
    vyco    =   zeros(Float64,(NC...)),
    vc      =   zeros(Float64,(NC...)),
    ΔTtop   =   zeros(Float64,NC.x),
    ΔTbot   =   zeros(Float64,NC.x),
)
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
    itmax   =   30000,              #   Maximum iterations
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
count       =   0
epsV0       =   0
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
    val     =   (E=zeros(NV.y),W=zeros(NV.y),S=zeros(NV.x),N=zeros(NV.x),
                    vyS=0.0,vyN=0.0,vxW=0.0,vxE=0.0),
)
# ----------------------------------------------------------------------- #
# Rayleigh Number ======================================================= #
# if P.Ra < 0
    P.Ra     =   P.ρ₀*P.g*P.α*P.ΔT*S.hsc^3/P.η₀/P.κ
# else
#     P.η₀     =   P.ρ₀*P.g*P.α*P.ΔT*S.hsc^3/P.Ra/P.κ
# end
filename    =   string("13_Blankenbach_",@sprintf("%.2e",P.Ra),
                        "_",NC.x,"_",NC.y,
                        "_",Ini.T)
# ----------------------------------------------------------------------- #
# Linear System of Equations ============================================ #
# Momentum Conservation Equation (MCE) ------
niterM =   50
ϵM     =   1e-8
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
KM      =   ExtendableSparseMatrix(ndof,ndof)      
δx      =   zeros(maximum(Num.Pt))
F       =   zeros(maximum(Num.Pt))
# Energy Conservation Equation (ECE) ------
niterT  =   10
ϵT      =   1e-12
K1      =   ExtendableSparseMatrix(ndof,ndof)
K2      =   ExtendableSparseMatrix(ndof,ndof)
rhs     =   zeros(ndof)
RT      =   zeros(NC...)
∂2T     =   (∂x2=zeros(NC.x, NC.y), ∂y2=zeros(NC.x, NC.y),
                ∂x20=zeros(NC.x, NC.y), ∂y20=zeros(NC.x, NC.y))
# ----------------------------------------------------------------------- #
# Time Loop ============================================================= #
for it = 1:T.itmax
    if it>1
        Time[it]  =   Time[it-1] + T.Δ
    end
    @printf("Time step: #%04d, Time [non-dim]: %04e\n ",it,
                    Time[it])
    # ------ MCE ------
    @printf("---Momentum Calculation ---\n")
    if it == 1
        D.vx[2:end-1,:]    .=  0.0
        D.vy[:,1:end-1]    .=  0.0
        D.Pt               .=  0.0
    end
    # Residual Calculation ------
    @. D.ρ  =   -P.Ra*D.T
    for iter = 1:niterM
        Residuals2Dc!(D,VBC,ε,τ,divV,Δ,1.0,1.0,Fm,FPt)
        F[Num.Vx]    =   Fm.x[:]
        F[Num.Vy]    =   Fm.y[:]
        F[Num.Pt]    =   FPt[:]
        @printf("||RM|| = %1.4e\n", norm(F)/length(F))
            norm(F)/length(F) < ϵM ? break : nothing
        # Update K ------
        KM      =   Assemblyc(NC, NV, Δ, 1.0, VBC, Num)
        # Solving Linear System of Equations ------
        δx      =   - KM \ F
        # Update unknown variables ------
        D.vx[:,2:end-1]     .+=  δx[Num.Vx]
        D.vy[2:end-1,:]     .+=  δx[Num.Vy]
        D.Pt                .+=  δx[Num.Pt]
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
    if mod(it,50) == 0 || it == T.itmax || it == 1
        p = heatmap(x.c,y.c,D.T',
                xlabel="x",ylabel="y",colorbar=true,
                title="Temperature",color=cgrad(:lajolla),
                aspect_ratio=:equal,xlims=(M.xmin, M.xmax),
                ylims=(M.ymin, M.ymax))
        contour!(p,x.c,y.c,D.T',lw=1,
                    color="white",cbar=false,alpha=0.5)
        quiver!(p,x.c2d[1:Pl.qinc:end,1:Pl.qinc:end],
            y.c2d[1:Pl.qinc:end,1:Pl.qinc:end],
            quiver=(D.vxc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc,
                    D.vyc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc),
            la=0.5,color="black",alpha=0.5)
        if save_fig == 1
            Plots.frame(anim)
        elseif save_fig == 0
            display(p)
        end
    end
    # ------------------------------------------------------------------- #
    # Advection ========================================================= #
    if it == 1
        @. D.vxco   =   D.vxc
        @. D.vyco   =   D.vyc
    end
    semilagc2D!(D.T,D.T_ex,D.vxc,D.vyc,D.vxco,D.vyco,x,y,T.Δ)
    @. D.vxco  =   D.vxc
    @. D.vyco  =   D.vyc
    # ------------------------------------------------------------------- #
    # Diffusion ========================================================= #
    @printf("---Energy Calculation---\n")
    # CNA2Dc!(D, 1.0, Δ.x, Δ.y, T.Δ, NC, TBC, rhs, K1, K2, Num)
    # Alternatively - defect correction ---
    # Update temperature field --- 
    @. D.T0     =   D.T
    for iter = 1:niter
        # Evaluate residual
        ComputeResiduals2Dc!( RT, D.T, D.T_ex, D.T0, D.T_ex0, ∂2T, 
                1.0, TBC, Δ, T.Δ[1]; C = 0.5 )
        @printf("||RT|| = %1.4e\n", norm(RT)/length(RT))
        norm(RT)/length(RT) < ϵT ? break : nothing
        # Assemble linear system
        K1  = AssembleMatrix2Dc(1.0, TBC, Num, NC, Δ, T.Δ[1];C=0.5)
        # Solve for temperature correction: Cholesky factorisation
        Kc = cholesky(K1.cscmatrix)
        # Solve for temperature correction: Back substitutions
        δT = -(Kc\RT[:])
        # Update temperature
        @. D.T += δT[Num.T]
    end
    # Update extended field for advection scheme --- 
    D.T_ex[2:end-1,2:end-1]     .=  D.T
    # ------------------------------------------------------------------- #
    # Nusselt Number ==================================================== #
    # Grid structure at the surface ---
    #   o - Centroids
    #   x - Vertices 
    #   □ - Ghost Nodes
    #
    #   □          □           □            □
    #   
    #        x --------- x --------- x
    #        |           |           |
    #   □    |     o     |     o     |      □ 
    #        |           |           |
    #        x --------- x --------- x      
    #        |           |           |
    #   □    |     o     |     o     |      □
    # --- 
    # Get temperature at the vertices 
    Tv1     =   zeros(NV.x,1)
    Tv2     =   zeros(NV.x,1)
    @. Tv1  =   (D.T_ex[1:end-1,end] + D.T_ex[2:end,end] + 
                    D.T_ex[1:end-1,end-1] + D.T_ex[2:end,end-1])/4
    @. Tv2  =   (D.T_ex[1:end-1,end-1] + D.T_ex[2:end,end-1] + 
                    D.T_ex[1:end-1,end-2] + D.T_ex[2:end,end-2])/4
    # Calculate temperature gradient --- 
    dTdy    =   zeros(NV.x,1)
    @. dTdy =   -(Tv1 - Tv2)/Δ.y
    # Calculate Nusselt number ---
    # Midpoint integration -
    # for i=1:NC.x
    #     Nus[it] +=   ((dTdy[i]+dTdy[i+1])/2)*Δ.x
    # end
    # Trapezoidal integration -
    for i = 1:NV.x
        if i == 1 || i == NV.x
            afac = 1
        else
            afac = 2
        end
        Nus[it]     += afac * dTdy[i]
    end
    Nus[it]     *=   Δ.x/2
    # Mean Temperature ---
    meanT[it,:] =   mean(D.T_ex,dims=1)
    # Root Mean Square Velocity ---
    meanV[it]   =   mean(D.vc)
    # ------------------------------------------------------------------- #
    # Check break ======================================================= #
    # If the maximum time is reached or if the models reaches steady 
    # state the time loop is stoped! 
    if Time[it] > 0.05
        epsC    =   0.001
        ind     =   findfirst(Time .> 
                        (Time[it] - 0.05))
        epsV    =   std(meanV[ind:it])
        if count == 0
            epsV0   =   epsV
            count   += 1
        end
        epsV        /= epsV0
        if save_fig == 1
            plot!(k,(it,log10((epsV))),
                xlabel="it",ylabel="log₁₀(εᵥ/eᵥ₀)",label="",
                markershape=:circle,markercolor=:black)
                    # display(k)
        end
        find    =   it
        @printf("ε_Vr = %g, ε_Cr = %g \n",epsV,epsC)
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
            xlabel="x",ylabel="y",colorbar=true,
            title="Temperature",color=cgrad(:lajolla),
            aspect_ratio=:equal,xlims=(M.xmin, M.xmax),
            ylims=(M.ymin, M.ymax))
    contour!(p2,x.c,y.c,D.T',lw=1,
                    color="white",cbar=false,
                    alpha=0.5)
    quiver!(p2,x.c2d[1:Pl.qinc:end,1:Pl.qinc:end],
            y.c2d[1:Pl.qinc:end,1:Pl.qinc:end],
            quiver=(D.vxc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc,
                    D.vyc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc),
            la=0.5,color="black",
            alpha=0.5)
if save_fig == 1
    savefig(k,string("./exercises/Correction/Results/13_BlankenbachBenchmark_iterations_",@sprintf("%.2e",P.Ra),
            "_",NC.x,"_",NC.y,
            "_",Ini.T,".png"))
    savefig(p2,string("./exercises/Correction/Results/13_BlankenbachBenchmark_Final_Stage_",@sprintf("%.2e",P.Ra),
            "_",NC.x,"_",NC.y,"_it_",find,"_",
            Ini.T,".png"))
elseif save_fig == 0
    display(p2)
    display(k)
end
# ----------------------------------------------------------------------- #
# Plot time serieses ==================================================== #
q2  =   plot(Time[1:find],Nus[1:find],
            xlabel="Time [ non-dim ]", ylabel="Nus",label="",
            layout=(2,1),suplot=1)
plot!(q2,Time[1:find],B.Nu[B.Nr].*ones(find,1),
            lw=0.5,color="red",linestyle=:dash,alpha=0.5,label="",
            layout=(2,1),suplot=1)
plot!(q2,Time[1:find],meanV[1:find],
            xlabel="Time [ non-dim ]", ylabel="V_{RMS}",label="",
            layout=(2,1),subplot=2)
plot!(q2,Time[1:find],B.Vrms[B.Nr].*ones(find,1),
            lw=0.5,color="red",linestyle=:dash,alpha=0.5,label="",
            layout=(2,1),subplot=2)
if save_fig == 1
    savefig(q2,string("./exercises/Correction/Results/13_BlankenbachBenchmark_TimeSeries_",@sprintf("%.2e",P.Ra),
                        "_",NC.x,"_",NC.y,"_",Ini.T,".png"))
elseif save_fig == 0
    display(q2)
end
# ----------------------------------------------------------------------- #
# Vertical temperature profiles ========================================= #
ind1    =   Int64(floor((M.xmax-M.xmin)/2/Δ.x)+1)
ind2    =   Int64(floor((M.xmax-M.xmin)/2/Δ.x))
Tprof   =   (D.T[ind1,:]+D.T[ind2,:])/2
q3  =   plot(Tprof,y.c,
        xlabel="T",ylabel="y",label=false)
scatter!(q3,(B.Tmin[B.Nr],-(1-B.ymin[B.Nr])),
            markershape=:rect,markersize=3,
            markercolor=:black,label=false)
scatter!(q3,(B.Tmax[B.Nr],-(1-B.ymax[B.Nr])),
            markershape=:rect,markersize=3,
            markercolor=:black,label=false)
if save_fig == 1
    savefig(q3,string("./exercises/Correction/Results/13_BlankenbachBenchmark_Profiles_",@sprintf("%.2e",P.Ra),
                        "_",NC.x,"_",NC.y,"_",Ini.T,".png"))
elseif save_fig == 0
    display(q3)
end
# ----------------------------------------------------------------------- #
end

BlankenbachBenchmark()
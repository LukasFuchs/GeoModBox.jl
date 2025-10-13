using Plots, ExtendableSparse
using GeoModBox
using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
using GeoModBox.AdvectionEquation.TwoD, GeoModBox.HeatEquation.TwoD
using GeoModBox.Scaling
using Statistics, LinearAlgebra
using Printf

function BlankenbachBenchmark(ncy,Ra,save_fig)
    # Define Initial Condition ========================================== #
    # Temperature - 
    #   1) circle, 2) gaussian, 3) block, 4) linear, 5) lineara
    # !!! Gaussian is not working!!! 
    Ini         =   (T=:lineara,)
    # ------------------------------------------------------------------- #
    # Plot Settings ===================================================== #
    Pl  =   (
        qinc        =   5,
        qsc         =   1.0e-3
    )
    # ------------------------------------------------------------------- #
    # Benchmark values ================================================== # 
    # Taken from Gerya (2019), Introduction to numerical geodynamic 
    if Ra == 1e4
        nr  =   1
    elseif Ra == 1e5
        nr = 2
    elseif Ra == 1e6
        nr = 3
    else 
        @printf("Error! Ra not defined for Benchmark")
        return
    end
    B       =   (
        Nr      =   nr,
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
    # ------------------------------------------------------------------- #
    # Geometry ========================================================== #
    M   =   Geometry(
        xmin    =   0.0,                #   [ m ] 
        xmax    =   1000e3,             #   [ m ]
        ymin    =   -1000e3,            #   [ m ]
        ymax    =   0.0,                #   [ m ]
    )
    # ------------------------------------------------------------------- #
    # Referenzparameter ================================================= #
    P   =   Physics(
        g       =   10.0,               #   Graviational Acceleration [m/s^2]
        ρ₀      =   4000.0,             #   Reference Density [kg/m^3]
        k       =   5.0,                #   Thermal Conductivity [ W/m/K ]
        cp      =   1250.0,             #   Heat capacity [ J/kg/K ]
        α       =   2.5e-5,             #   THermal Expansion [ K^-1 ]
        ΔT      =   1000.0,             #   Temperature Difference
        # If Ra < 0, Ra will be calculated from the reference parameters.
        # If Ra is defined, the reference viscosity will be adjusted such that
        # the scaling parameters result in the given Ra.
        Ra      =   Ra,              #   Rayleigh number
    )
    # ------------------------------------------------------------------- #
    # Define Scaling Constants ========================================== # 
    S   =   ScalingConstants!(M,P)
    # ------------------------------------------------------------------- #
    # Numerical Grid ==================================================== #
    NC  =   (
        x   =   ncy,
        y   =   ncy,
    )
    NV      =   (
        x   =   NC.x + 1,
        y   =   NC.y + 1,
    )
    Δ       =   GridSpacing(
        x   =   (M.xmax - M.xmin)/NC.x,
        y   =   (M.ymax - M.ymin)/NC.y,
    )
    # ------------------------------------------------------------------- #
    # Data Arrays ======================================================= #
    D       =   DataFields(
        Q       =   zeros(Float64,(NC...)),
        T       =   zeros(Float64,(NC...)),
        T0      =   zeros(Float64,(NC...)),
        T_ex    =   zeros(Float64,(NC.x+2,NC.y+2)),
        ρ       =   ones(Float64,(NC...)),
        vx      =   zeros(Float64,(NV.x,NV.y+1)),
        vy      =   zeros(Float64,(NV.x+1,NV.y)),    
        Pt      =   zeros(Float64,(NC...)),
        vxc     =   zeros(Float64,(NC...)),
        vyc     =   zeros(Float64,(NC...)),
        vc      =   zeros(Float64,(NC...)),
        ΔTtop   =   zeros(Float64,NC.x),
        ΔTbot   =   zeros(Float64,NC.x),
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
    # Residuals ---
    Fm     =    (
        x       =   zeros(Float64,NV.x, NC.y), 
        y       =   zeros(Float64,NC.x, NV.y)
    )
    FPt         =   zeros(Float64,NC...)
    # ------------------------------------------------------------------- #
    # Time parameters =================================================== #
    T   =   TimeParameter(
        tmax    =   1000000.0,          #   [ Ma ]
        Δfacc   =   0.9,                #   Courant time factor
        Δfacd   =   0.9,                #   Diffusion time factor
        itmax   =   10000,              #   Maximum iterations
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
    # ------------------------------------------------------------------- #
    # Scaling laws ====================================================== #
    ScaleParameters!(S,M,Δ,T,P,D)
    # ------------------------------------------------------------------- #
    # Coordinates ======================================================= #
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
    # ------------------------------------------------------------------- #
    # Initial Condition ================================================= #
    # Temperature ------
    IniTemperature!(Ini.T,M,NC,D,x,y;Tb=P.Tbot,Ta=P.Ttop)
    # ------------------------------------------------------------------- #
    # Boundary Conditions =============================================== #
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
    # ------------------------------------------------------------------- #
    # Reference Viscosity =============================================== #
    P.η₀     =   P.ρ₀*P.g*P.α*P.ΔT*S.hsc^3/P.Ra/P.κ
    # ------------------------------------------------------------------- #
    # Linear System of Equations ======================================== #
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
    # ------------------------------------------------------------------- #
    # Time Loop ========================================================= #
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
        # Calculate time stepping ======================================= #
        T.Δc        =   T.Δfacc * minimum((Δ.x,Δ.y)) / 
                (sqrt(maximum(abs.(D.vx))^2 + maximum(abs.(D.vy))^2))
        T.Δd        =   T.Δfacd * (1.0 / (2.0 *(1.0/Δ.x^2 + 1/Δ.y^2)))
        T.Δ         =   minimum([T.Δd,T.Δc])
        if Time[it] > T.tmax
            T.Δ         =   T.tmax - Time[it-1]
            Time[it]    =   Time[it-1] + T.Δ
            it          =   T.itmax
        end
        # Advection ===================================================== #
        semilagc2D!(D.T,D.T_ex,D.vxc,D.vyc,[],[],x,y,T.Δ)
        # --------------------------------------------------------------- #
        # Diffusion ===================================================== #
        CNA2Dc!(D, 1.0, Δ.x, Δ.y, T.Δ, D.ρ, 1.0, NC, TBC, rhs, K1, K2, Num)
        # --------------------------------------------------------------- #
        # Nusselt Number ================================================ #
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
        # --------------------------------------------------------------- #
        # Check break =================================================== #
        # If the maximum time is reached or if the models reaches steady 
        # state the time loop is stoped! 
        if Time[it] > 0.0038
            epsC    =   1e-3; 
            ind     =   findfirst(Time .> 
                            (Time[it] - 0.0038))
            epsV    =   std(meanV[ind:it])
            find    =   it
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
    # Save final figure ================================================= #
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
        savefig(p2,string("./exercises/Correction/Results/13_BlankenbachBenchmark_ResTest_Final_Stage_",@sprintf("%.2e",P.Ra),
                "_",NC.x,"_",NC.y,"_it_",find,"_",
                Ini.T,".png"))
    elseif save_fig == 0
        display(p2)
    end
    # ------------------------------------------------------------------- #
    # Plot time serieses ================================================ #
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
        savefig(q2,string("./exercises/Correction/Results/13_BlankenbachBenchmark_RestTest_TimeSeries_",@sprintf("%.2e",P.Ra),
                            "_",NC.x,"_",NC.y,"_",Ini.T,".png"))
    elseif save_fig == 0
        display(q2)
    end
    # ------------------------------------------------------------------- #
    return Nus[find], meanT[find], meanV[find]
end

# ======================================================================= #
# Run Resolution Test =================================================== #
# ======================================================================= #
start=time()
save_fig    =   1
# Rayleigh Number ======================================================= #
#   Here, only for 1e4, 1e5, 1e6
Ra      =   1e4
# ----------------------------------------------------------------------- #
# Benchmark Values ====================================================== # 
# Taken from Gerya (2019), Introduction to numerical geodynamic 
if Ra == 1e4
    nr  =   1
elseif Ra == 1e5
    nr = 2
elseif Ra == 1e6
    nr = 3
else 
    @printf("Error! Ra not defined for Benchmark")
    return
end
B       =   (
    # Nr      =   nr,
    # Nusselt Number at the top
    Nu      =   [4.8844,10.534,21.972,10.066,6.9229],   
    # Root Mean Square Velocity-scaled
    Vrms    =   [42.865,193.21,833.99,480.43,171.76],   
    # Non-dimensional temper^ature gradients in the model corner
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
# Number of grid points in upper thermal boundary layer ================= #
n   =   [2,3,4,5,6,7,8,9,10]
# ----------------------------------------------------------------------- #
# Grid Resultion ======================================================== #
ncy     =   zeros(Int64,length(n))
@. ncy  =   round((n-1)*(Ra/4)^(1/3))
# @show ncy
# ----------------------------------------------------------------------- #
# Statistical values ==================================================== #
Nus     =   zeros(length(n))
meanT   =   zeros(length(n))
meanV   =   zeros(length(n))
# Resolution Loop ======================================================= #
for k in eachindex(ncy)
    Nus[k],meanT[k],meanV[k] = BlankenbachBenchmark(ncy[k],Ra,save_fig)
end
# @show Nus, meanT, meanV
# ----------------------------------------------------------------------- #
xmin    =   1e-5
xmax    =   1e-1
p1  =   scatter((1 ./ncy./ncy,Nus),
            markershape=:circle,markersize=4,
            markercolor=:black,label="",
            xlims=[xmin,xmax],ylims=[B.Nu[nr]*0.5,B.Nu[nr]*1.3],
            xscale=:log10,yscale=:log10,
            xlabel="1/nx/ny",ylabel="Nu",
            layout=(3,1),subplot=1)
plot!(p1,LinRange(xmin,1,10),B.Nu[nr].*ones(10,1),
            lw=0.5,color="red",linestyle=:dash,alpha=0.5,label="",
            layout=(3,1),suplot=1)
scatter!(p1,(1 ./ncy./ncy,meanV),
            markershape=:circle,markersize=4,
            markercolor=:black,label="",
            xlims=[xmin,xmax],ylims=[B.Vrms[nr]*0.5,B.Vrms[nr]*2.0],
            xscale=:log10,yscale=:log10,
            xlabel="1/nx/ny",ylabel="V_{RMS}",
            layout=(3,1),subplot=2)
plot!(p1,LinRange(xmin,1,10),B.Vrms[nr].*ones(10,1),
            lw=0.5,color="red",linestyle=:dash,alpha=0.5,label="",
            layout=(3,1),subplot=2)
scatter!(p1,(1 ./ncy./ncy,meanT),
            markershape=:circle,markersize=4,
            markercolor=:black,label="",
            xlims=[xmin,xmax],ylims=[1e-1,10],
            xscale=:log10,yscale=:log10,
            xlabel="1/nx/ny",ylabel="⟨T⟩",
            layout=(3,1),subplot=3)

if save_fig == 1
    savefig(p1,string("./exercises/Correction/Results/13_BlankenbachBenchmark_RestTest_Ra_",
                        @sprintf("%.2e",Ra),".png"))
elseif save_fig == 0
    display(p1)
end
stop=time()
@printf("Runtime: %g s",stop-start)

using Plots
using ExtendableSparse
using GeoModBox
using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
using GeoModBox.AdvectionEquation.TwoD, GeoModBox.Scaling
using GeoModBox.Tracers.TwoD
using Base.Threads
using Printf, LinearAlgebra
using TimerOutputs
using Statistics

function VanKekenBenchmark()

    to          =   TimerOutput()
    @timeit to "Ini" begin
    save_fig    =   1
    # Define Initial Condition ========================================== #
    Ini         =   (p=:RTI,) 
    avg         =   :arith
    # ------------------------------------------------------------------- #
    # Plot Settings ===================================================== #
    Pl  =   (
        qinc    =   10,
        mainc   =   5,
        qsc     =   0.1, # 100*(60*60*24*365.25)*3
    )
    # ------------------------------------------------------------------- #
    # Geometry ========================================================== #
    M       =   Geometry(
        ymin    =   -4.0e3,     # [ m ]
        ymax    =   0.0,
        xmin    =   0.0,
    )
    # Dimensionless benchmark ratios
    λratio      =   2.0 * 0.9142
    Aratio      =   0.02
    hc          =   0.8

    Hdim        =   M.ymax - M.ymin
    M.xmax      =   0.5 * λratio * Hdim
    # -------------------------------------------------------------------- #
    # Physics ============================================================ #
    P   =   Physics(
        g       =   10.0,                #   Gravitational acceleration [ m/s^2 ]
        κ       =   1e-6,
        η₀      =   1e20,                #   Viscosity composition 0 [ Pa s ]
        ρ₀      =   3000.0,              #   Density composition 0 [ kg/m^3 ]
    )
    # 0 - upper layer; 1 - lower layer
    # η₀      =   1e20                #   Viscosity composition 0 [ Pa s ]
    η₁      =   1e20                #   Viscosity composition 1 [ Pa s]
    ηᵣ      =   log10(η₁/P.η₀)
    η       =   [P.η₀,η₁]             #   Viscosity for phases 

    h₀      =   hc * (M.ymax-M.ymin)        # Thickness upper layer
    h₁      =   (1-hc) * (M.ymax-M.ymin)    # Thickness lower layer
    yc      =   abs(M.ymin + h₁)            # Layer interface depth

    Rb      =   1

    phase   =   [0,1]
    # -------------------------------------------------------------------- #
    # Define Scaling Constants ============================================== # 
    S   =   ScalingConstants!(M,P)
    # ----------------------------------------------------------------------- #
    # Grid =============================================================== # 
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
    # -------------------------------------------------------------------- #
    # Animation and Plot Settings ======================================= #
    path        =   string("./examples/StokesEquation/2D/Results/")
    anim        =   Plots.Animation(path, String[] )
    filename    =   string("VanKeKen_Benchmark_ηr_",round(ηᵣ),
                            "_tracers_DC_scaled_",avg)
    # ------------------------------------------------------------------- #
    # Allocation ======================================================== #
    D       =   DataFields(
        # Q       =   zeros(Float64,(NC...)),
        p       =   zeros(Float64,(NC...)),
        p_ex    =   zeros(Float64,(NC.x+2,NC.y+2)),
        ρ       =   zeros(Float64,(NC...)),  
        ρ_ex    =   zeros(Float64,(NC.x+2,NC.y+2)),    
        cp      =   zeros(Float64,(NC...)),
        vx      =   zeros(Float64,(NV.x,NV.y+1)),
        vy      =   zeros(Float64,(NV.x+1,NV.y)),    
        Pt      =   zeros(Float64,(NC...)),
        vxc     =   zeros(Float64,(NC...)),
        vyc     =   zeros(Float64,(NC...)),
        vc      =   zeros(Float64,(NC...)),
        wt      =   zeros(Float64,(NC.x,NC.y)),
        wte     =   zeros(Float64,(NC.x+2,NC.y+2)),
        wtv     =   zeros(Float64,(NV.x,NV.y)),
        ηc      =   zeros(Float64,NC...),
        η_ex    =   zeros(Float64,(NC.x+2,NC.y+2)),
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
    # ------------------------------------------------------------------- #
    # Boundary Conditions =============================================== #
    VBC     =   (
        type    =   (E=:freeslip,W=:freeslip,S=:noslip,N=:noslip),
        val     =   (E=zeros(NV.y),W=zeros(NV.y),S=zeros(NV.x),N=zeros(NV.x),
                        vyS=0.0,vyN=0.0,vxW=0.0,vxE=0.0),
    )
    # ------------------------------------------------------------------- #
    # Time ============================================================== #
    T   =   TimeParameter(
        tmax    =   5000.0,         #   [ Ma ]
        Δfacc   =   1.0,            #   Courant time factor
        itmax   =   500,            #   Maximum iterations; 50
    )
    T.tmax      =   T.tmax*1e6*T.year    #   [ s ]
    T.Δ         =   0.0

    Time        =   zeros(T.itmax)
    V_RMS       =   zeros(T.itmax)
    final_step  =   0
    # ------------------------------------------------------------------- #
    # Scaling laws ========================================================== #
    ScaleParameters!(S,M,Δ,T,P,D)
    # Correct parameters
    # Scaled geometrical parameters
    H  = M.ymax - M.ymin

    h₀ = hc * H
    h₁ = (1.0 - hc) * H

    yc = abs(M.ymin + h₁)
    λ  = λratio * H
    δA = - Aratio * H

    η   ./=     P.η₀
    # ----------------------------------------------------------------------- #
    # Coordinates ======================================================= #
    x       =   (
        c   =   LinRange(M.xmin+Δ.x/2,M.xmax-Δ.x/2,NC.x),
        ce  =   LinRange(M.xmin - Δ.x/2.0, M.xmax + Δ.x/2.0, NC.x+2),
        v   =   LinRange(M.xmin,M.xmax,NV.x),
    )
    y       =   (
        c   =   LinRange(M.ymin+Δ.y/2,M.ymax-Δ.y/2,NC.y),
        ce  =   LinRange(M.ymin - Δ.y/2.0, M.ymax + Δ.y/2.0, NC.y+2),
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
    # Tracer Advection ================================================== #
    @timeit to "Tracer Ini" begin
    nmx,nmy     =   5,5
    noise       =   1
    nmark       =   nmx*nmy*NC.x*NC.y
    Aparam      =   :phase
    MPC         =   (
            c       =   zeros(Float64,(NC.x,NC.y)),
            v       =   zeros(Float64,(NV.x,NV.y)),
            th      =   zeros(Float64,(nthreads(),NC.x,NC.y)),
            thv     =   zeros(Float64,(nthreads(),NV.x,NV.y)),
    )
    MAVG        = (
            PC_th   =   [similar(D.wte) for _ = 1:nthreads()],  # per thread
            PV_th   =   [similar(D.ηv) for _ = 1:nthreads()],   # per thread
            wte_th  =   [similar(D.wte) for _ = 1:nthreads()],  # per thread
            wtv_th  =   [similar(D.wtv) for _ = 1:nthreads()],  # per thread
    )
    Ma      =   IniTracer2D(Aparam,nmx,nmy,Δ,M,NC,noise,Ini.p,phase;
                            xc=0.0,yc=yc,λ,δA)
    # RK4 weights ---
    rkw     =   1.0/6.0*[1.0 2.0 2.0 1.0]   # for averaging
    rkv     =   1.0/2.0*[1.0 1.0 2.0 2.0]   # for time stepping
    # Count marker per cell ---
    CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,NV)
    # Interpolate from markers to cell ---
    Markers2Cells(Ma,nmark,MAVG.PC_th,D.p_ex,MAVG.wte_th,D.wte,x,y,Δ,:phase,phase)
    D.p .= D.p_ex[2:end-1, 2:end-1]
    # Phase field (0 or 1) -> nondimensional density anomaly
    @. D.ρ = -Rb * D.p
    Markers2Cells(Ma,nmark,MAVG.PC_th,D.η_ex,MAVG.wte_th,D.wte,x,y,Δ,Aparam,η;avgm=avg)
    D.ηc    .=   D.η_ex[2:end-1,2:end-1]
    Markers2Vertices(Ma,nmark,MAVG.PV_th,D.ηv,MAVG.wtv_th,D.wtv,x,y,Δ,Aparam,η;avgm=avg)
    end
    # ------------------------------------------------------------------- #
    # System of Equations =============================================== #
    # Iterations
    niter       =   50
    atol        =   1e-8        #   Absolute tolerance
    rtol        =   1e-5        #   Relative tolerance; r = atolM0/atolM
    RM          =   0.0         #   Initialize absolute residual    
    RMrel       =   0.0         #   Initialize relative residual 
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
    # Residuals ---
    Fm     =    (
        x       =   zeros(Float64,NV.x, NC.y), 
        y       =   zeros(Float64,NC.x, NV.y)
    )
    FPt     =   zeros(Float64,NC...)      
    # ------------------------------------------------------------------- #
    end
    # Time Loop ========================================================= #
    @timeit to "Time Loop" begin
    for it = 1:T.itmax
        R0      =   0.0
        # Update Time ---
        if it > 1
            Time[it]   =   Time[it-1] + T.Δ 
        end
        @printf("Time step: #%04d, Time: %04e\n ",it,
                    Time[it])
        # Momentum Equation ===
        # Initial Residual ---------------------------------------------- #
        if it == 1
            D.vx    .=  0.0
            D.vy    .=  0.0
            D.Pt    .=  0.0
        end
        @timeit to "Solution Iteration" begin
        # Assemble Coefficients ========================================= #
        @timeit to "Assembly" begin
        K       =   Assembly(NC, NV, Δ, D.ηc, D.ηv, VBC, Num)
        Kfac    =   lu(K.cscmatrix)
        end
        for iter=1:niter
            @timeit to "Residual" begin
            Residuals2D!(D,VBC,ε,τ,divV,Δ,D.ηc,D.ηv,1.0,Fm,FPt)
            F[Num.Vx]   .=   Fm.x
            F[Num.Vy]   .=   Fm.y
            F[Num.Pt]   .=   FPt
            RM          =   norm(F)/length(F)
            if iter == 1
                R0 = max(RM, eps(Float64))
            end
            RMrel       =   RM/R0
            @printf("   MCE %2d: ||R|| = %1.4e, ||R||/||R₀|| = %1.4e\n",iter,RM,RMrel)
            (RM < atol || RM/R0 < rtol) && break
            end
            # --------------------------------------------------------------- #
            # Solution of the linear system ================================= #
            # @printf("hello")
            @timeit to "Solution" begin
            δx      =   - (Kfac \ F)
            end
            # --------------------------------------------------------------- #
            # Update Unknown Variables ====================================== #
            D.vx[:,2:end-1]     .+=  δx[Num.Vx]
            D.vy[2:end-1,:]     .+=  δx[Num.Vy]
            D.Pt                .+=  δx[Num.Pt]
        end
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
        # ---
        if mod(it,4) == 0 || final_step == 1 || it == 1
            p = heatmap(x.c,y.c,D.p',color=:inferno,
                        xlabel="x",ylabel="y",colorbar=true,
                        title="Phase",
                        aspect_ratio=:equal,xlims=(M.xmin, M.xmax),                             
                        ylims=(M.ymin, M.ymax),
                        layout=(2,2),subplot=1)
            scatter!(p,Ma.x[1:Pl.mainc:end],Ma.y[1:Pl.mainc:end],
                        ms=1,ma=0.5,mc=Ma.phase[1:Pl.mainc:end],markerstrokewidth=0.0,
                        xlabel="x",ylabel="y",colorbar=true,
                        title="tracers",label="",
                        aspect_ratio=:equal,xlims=(M.xmin, M.xmax), 
                        ylims=(M.ymin, M.ymax),
                        layout=(2,2),subplot=2)
            heatmap!(p,x.c,y.c,D.vc',
                        xlabel="x",ylabel="y",colorbar=true,
                        title="V_c",color=cgrad(:batlow),
                        aspect_ratio=:equal,xlims=(M.xmin, M.xmax),
                        ylims=(M.ymin, M.ymax),
                        layout=(2,2),subplot=4)
            quiver!(p,x.c2d[1:Pl.qinc:end,1:Pl.qinc:end],
                        y.c2d[1:Pl.qinc:end,1:Pl.qinc:end],
                        quiver=(D.vxc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc,
                                D.vyc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc),        
                        la=0.5,color="white",layout=(2,2),subplot=4)
            heatmap!(p,x.c,y.c,log10.(D.ηc'),color=reverse(cgrad(:roma)),
                        xlabel="x",ylabel="y",title="log10(η_c)",
                        # clims=(15,27),
                        aspect_ratio=:equal,xlims=(M.xmin, M.xmax), 
                        ylims=(M.ymin, M.ymax),colorbar=true,
                        layout=(2,2),subplot=3)
            if save_fig == 1
                Plots.frame(anim)
            elseif save_fig == 0
                display(p)
            end
        end
        if final_step == 1
            @printf(" Maximum Time reached!\n")
            break
        end
        # Calculate Time Stepping ---
        T.Δ        =   T.Δfacc * minimum((Δ.x,Δ.y)) / 
                (sqrt(maximum(abs.(D.vx))^2 + maximum(abs.(D.vy))^2))
        if Time[it] >= T.tmax
            T.Δ         =   T.tmax - Time[it-1]
            Time[it]    =   Time[it-1] + T.Δ
            final_step  =   1
        end
        # Advection ===
        @timeit to "Tracer Advection" begin
        # Advect tracers ---
        @printf("Running on %d thread(s)\n", nthreads())  
        AdvectTracer2D(Ma,nmark,D,x,y,T.Δ,Δ,NC,rkw,rkv;style=1)
        CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,NV)
        # Interpolate phase from tracers to grid ---
        Markers2Cells(Ma,nmark,MAVG.PC_th,D.p_ex,MAVG.wte_th,D.wte,x,y,Δ,Aparam,phase)
        D.p     .=  D.p_ex[2:end-1,2:end-1]
        # Nondimensional density anomaly
        # Phase field (0 or 1) -> nondimensional density anomaly
        @. D.ρ = -Rb * D.p
        Markers2Cells(Ma,nmark,MAVG.PC_th,D.η_ex,MAVG.wte_th,D.wte,x,y,Δ,Aparam,η;avgm=avg)
        D.ηc    .=   D.η_ex[2:end-1,2:end-1]
        Markers2Vertices(Ma,nmark,MAVG.PV_th,D.ηv,MAVG.wtv_th,D.wtv,x,y,Δ,Aparam,η;avgm=avg)
        end
        V_RMS[it]   =   sqrt(mean(D.vxc.^2 .+ D.vyc.^2))
    end # End Time Loop
    end
    # Save Animation ---
    if save_fig == 1
        # Write the frames to a GIF file
        Plots.gif(anim, string( path, filename, ".gif" ), fps = 15)
        foreach(rm, filter(startswith(string(path,"00")), readdir(path,join=true)))
    end
    # Plot time serieses ==================================================== #
    p2  =   plot(Time,V_RMS,
                xlabel="Time", ylabel="V_RMS",label="")
    if save_fig == 1
        savefig(p2,string("./examples/StokesEquation/2D/Results/VanKekenBenchmark_TimeSeries_Scaled",
                            "_",NC.x,"_",NC.y,"_",avg,".png"))
    elseif save_fig == 0
        display(p2)
    end
    display(to)
end

VanKekenBenchmark()
using Plots
using ExtendableSparse
using GeoModBox
using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
using GeoModBox.AdvectionEquation.TwoD
using GeoModBox.Tracers.TwoD
using Base.Threads
using Printf, LinearAlgebra

function RTI_GrowthRate()
    plot_fields     =:no
    save_fig        = 1
    Pl  =   (
        qinc    =   5, 
        qsc     =   100*(60*60*24*365.25)*5e1,
    )
    # Define Initial Condition ========================================== #
    addnoise    =   [0 1]
    # Density Averaging ---
    #   centroids or vertices
    ρavg        =   :centroids
    nm          =   [2 4 6 8 10 12 14 16 18 20]
    # nc          =   
    # Initial Marker distribution ---
    Ini         =   (p=:RTI,) 
    # Perturbation wavelength [ m ]
    λᵣ          =   4e3
    # ------------------------------------------------------------------- #
    # Physics =========================================================== #
    g           =   9.81                    #   Gravitational acceleration [ m/s^2 ]
    # Lower layer ---
    ρ₁          =   2900.0                  #   Density composition 1 [ kg/m^3 ]
    η₁          =   1e19                    #   Viscosity composition 1 [ Pa s]
    # Upper layer --- 
    ρ₀          =   3000.0                  #   Density composition 0 [ kg/m^3 ]
    # ---
    ρ           =   [ρ₀,ρ₁]                 #   Density for phases
    ηᵣ          =   [1e-6 1 500]            #   Viscosity ratio
    phase       =   [0,1]
    # ------------------------------------------------------------------- #
    # Divisional factor of the amplitude following Gerya (2009) --------- #
    delfac      =   [15 150 1500]
    ms          =   zeros(3)
    ms          =   [3]
    mc          =   ["black","red","yellow"] # δA
    # Plot Settings ===================================================== #
    q   =   plot(layout=(size(addnoise,2),size(ηᵣ,2)))
    # ------------------------------------------------------------------- #
    # Geometry ========================================================== #
    M       =   Geometry(
        ymin    =   -3.0e3,     #   [ m ]
        ymax    =   0.0,
        xmin    =   0.0,
    )
    # Perturbation properties ---
    λ           =   λᵣ[1]       #   [ m ]
    # ---
    ar          =   Int64(round(2 * λ / (M.ymax-M.ymin)))  # aspect ratio
    M.xmax      =   (M.ymax-M.ymin)*ar
    @printf("   xmax: %g \n",M.xmax)
    @printf("   λ = %g\n",λ)
    # ------------------------------------------------------------------- #
    # Grid ============================================================== # 
    NC  =   (
        x   =   50,
        y   =   50,
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
    # Allocation ======================================================== #
    D       =   (
        ρ       =   zeros(Float64,(NC...)),
        ρv      =   zeros(Float64,NV...),
        ρe      =   zeros(Float64,(NC.x+2,NC.y+2)),
        p       =   zeros(Float64,(NC...)),
        cp      =   zeros(Float64,(NC...)),
        vx      =   zeros(Float64,(NV.x,NV.y+1)),
        vy      =   zeros(Float64,(NV.x+1,NV.y)),    
        Pt      =   zeros(Float64,(NC...)),
        vxc     =   zeros(Float64,(NC...)),
        vyc     =   zeros(Float64,(NC...)),
        vc      =   zeros(Float64,(NC...)),
        wte     =   zeros(Float64,(NC.x+2,NC.y+2)),
        wtv     =   zeros(Float64,(NV.x,NV.y)),
        ηc      =   zeros(Float64,NC...),
        ηce     =   zeros(Float64,(NC.x+2,NC.y+2)),
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
        val     =   (E=zeros(NV.y),W=zeros(NV.y),S=zeros(NV.x),N=zeros(NV.x)),
    )
    # Analytical Solution =============================================== #
    λₐ          =   λᵣ                      # [ m ]
    ϕ₁          =   0.0
    ϕ₂          =   0.0
    c11         =   0.0
    d12         =   0.0
    i21         =   0.0
    j22         =   0.0
    # Arrays ---
    PP  =   (
        ϕ       =   [0.0],
        K       =   [0.0],
        Q       =   [0.0],
        ϕₐ      =   [0.0],
        Kₐ      =   zeros(1,length(ηᵣ)),
        ε       =   [0.0],
    )
    ϕ₁          =   (2*π*((M.ymax-M.ymin)/2))/λₐ
    ϕ₂          =   (2*π*((M.ymax-M.ymin)/2))/λₐ
    # ------------------------------------------------------------------- #
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
    # System of Equations =============================================== #
    # Iterations
    niter       =   10
    ϵ           =   1e-10
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
    for n in eachindex(addnoise)
        for o in eachindex(ηᵣ) # Loop over viscosity ratio
            @printf("    ηᵣ = %g\n",ηᵣ[o])
            # Physics =================================================== #
            # 0 - upper layer; 1 - lower layer
            η₀      =   η₁*ηᵣ[o]    #   Viscosity composition 0 [ Pa s ]
            η       =   [η₀,η₁]     #   Viscosity for phases 
            @printf("    η₀ = %g\n",η₀)
            # ----------------------------------------------------------- #
            # Analytical Solution ======================================= #
            c11     =   (η₀*2*ϕ₁^2)/
                            (η₁*(cosh(2*ϕ₁) - 1 - 2*ϕ₁^2)) - 
                            (2*ϕ₂^2)/
                            (cosh(2*ϕ₂) - 1 - 2*ϕ₂^2)
            d12     =   (η₀*(sinh(2*ϕ₁) - 2*ϕ₁))/
                            (η₁*(cosh(2*ϕ₁) - 1 - 2*ϕ₁^2)) + 
                            (sinh(2*ϕ₂) - 2*ϕ₂)/
                            (cosh(2*ϕ₂) - 1 - 2*ϕ₂^2)
            i21     =   (η₀*ϕ₂*(sinh(2*ϕ₁) + 2*ϕ₁))/
                            (η₁*(cosh(2*ϕ₁) - 1 - 2*ϕ₁^2)) + 
                            (ϕ₂*(sinh(2*ϕ₂) + 2*ϕ₂))/
                            (cosh(2*ϕ₂) - 1 - 2*ϕ₂^2)
            j22     =   (η₀*2*ϕ₁^2*ϕ₂)/
                            (η₁*(cosh(2*ϕ₁) - 1 - 2*ϕ₁^2)) - 
                            (2*ϕ₂^3)/
                            (cosh(2*ϕ₂) - 1 - 2*ϕ₂^2)
            
            PP.Kₐ[o]    =   -d12/(c11*j22 - d12*i21)
            PP.ϕₐ[1]    =   ϕ₁
            # ----------------------------------------------------------- #
            for l in eachindex(delfac) # Loop over perturbation amplitude
                δA          =   -(M.ymax-M.ymin)/2/delfac[l]    #   Amplitude [ m ]
                @printf("    δA = %g\n",δA)
                for k in eachindex(nm) # Loop over marker numbers
                    @printf("    nm = %g\n",nm[k])
                    # Tracer Advection ================================== #
                    nmx,nmy =   nm[k],nm[k]
                    noise   =   addnoise[n]
                    nmark   =   nmx*nmy*NC.x*NC.y
                    Aparam  =   :phase
                    # --------------------------------------------------- #
                    # Initialize Tracer Position ---
                    Ma      =   IniTracer2D(Aparam,nmx,nmy,Δ,M,NC,noise,Ini.p,phase;λ,δA)
                    # Count tracer per cell ---
                    CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,NV,1)
                    # Interpolate density --- 
                    if ρavg==:centroids
                        # Interpolate density from markers to cell ---
                        Markers2Cells(Ma,nmark,MAVG.PC_th,D.ρe,MAVG.wte_th,D.wte,x,y,Δ,Aparam,ρ)
                        D.ρ     .=   D.ρe[2:end-1,2:end-1]  
                    elseif ρavg==:vertices 
                        # Interpolate density from markers to vertices ---
                        Markers2Vertices(Ma,nmark,MAVG.PV_th,D.ρv,MAVG.wtv_th,D.wtv,x,y,Δ,Aparam,ρ)
                        for i = 1:NC.x
                            D.ρ[i,:]    .=   (D.ρv[i,1:end-1] .+ 
                                                D.ρv[i,2:end,:] .+ 
                                                D.ρv[i+1,1:end-1] .+ 
                                                D.ρv[i+1,2:end])./4
                        end
                    end     
                    # Interpolate Viscosity ---
                    Markers2Cells(Ma,nmark,MAVG.PC_th,D.ηce,MAVG.wte_th,D.wte,x,y,Δ,Aparam,η)
                    D.ηc    .=   D.ηce[2:end-1,2:end-1]
                    Markers2Vertices(Ma,nmark,MAVG.PV_th,D.ηv,MAVG.wtv_th,D.wtv,x,y,Δ,Aparam,η)
                    # --------------------------------------------------- #
                    # ------------------------------------------------------- #
                    # Momentum Equation ===
                    D.vx    .=  0.0
                    D.vy    .=  0.0
                    D.Pt    .=  0.0
                    @. δx   =   0.0
                    @. F    =   0.0
                    for iter=1:niter
                        # Initial Residual -------------------------------------- #
                        Residuals2D!(D,VBC,ε,τ,divV,Δ,D.ηc,D.ηv,g,Fm,FPt)
                        F[Num.Vx]   =   Fm.x[:]
                        F[Num.Vy]   =   Fm.y[:]
                        F[Num.Pt]   =   FPt[:]
                        @printf("||R|| = %1.4e\n", norm(F)/length(F))
                        norm(F)/length(F) < ϵ ? break : nothing
                        # Assemble Coefficients ================================= #
                        K       =   Assembly(NC, NV, Δ, D.ηc, D.ηv, VBC, Num)
                        # ------------------------------------------------------- #
                        # Solution of the linear system ========================= #
                        δx      =   - K \ F
                        # ------------------------------------------------------- #
                        # Update Unknown Variables ============================== #
                        D.vx[:,2:end-1]     .+=  δx[Num.Vx]
                        D.vy[2:end-1,:]     .+=  δx[Num.Vy]
                        D.Pt                .+=  δx[Num.Pt]
                    end
                    # ------------------------------------------------------- #
                    # Get the velocity on the centroids ---
                    # Just for visualization purposes
                    for i = 1:NC.x
                        for j = 1:NC.y
                            D.vxc[i,j]  = (D.vx[i,j+1] + D.vx[i+1,j+1])/2
                            D.vyc[i,j]  = (D.vy[i+1,j] + D.vy[i+1,j+1])/2
                        end
                    end
                    @. D.vc        = sqrt(D.vxc^2 + D.vyc^2)
                    # ---
                    # Calculate diapir growth rate ---
                    xwave       =   (M.xmax-M.xmin)/2  
                    ywave       =   (M.ymax-M.ymin)/2 + δA

                    xn          =   Int64(floor((xwave+Δ.x/2)/Δ.x))
                    yn          =   Int64(floor(((M.ymax-M.ymin)-ywave)/Δ.y)) + 1

                    dx          =   (xwave+Δ.x/2)/Δ.x - xn
                    dy          =   abs(((M.ymax-M.ymin)-ywave)/Δ.y - yn)

                    wvy     =   (1.0-dx)*(1.0-dy) * D.vy[xn+1,yn] + 
                                    dx*(1.0-dy) * D.vy[xn+2,yn] + 
                                    (1.0-dx)*dy * D.vy[xn+1,yn+1] + 
                                    dx*dy * D.vy[xn+2,yn+1]

                    PP.Q[1] =   (ρ₀-ρ₁)*(M.ymax-M.ymin)/2.0*g/2.0/η₁
                    PP.K[1] =   abs(wvy)/abs(δA)/PP.Q[1]
                    PP.ϕ[1] =   2*π*(M.ymax-M.ymin)/2/λ
                    PP.ε[1] =   abs((PP.Kₐ[o]-PP.K[1])/PP.Kₐ[o])*100.0
                    
                    if plot_fields==:yes
                        p = heatmap(x.c./1e3,y.c./1e3,D.ρ',color=:inferno,
                                    xlabel="x[km]",ylabel="y[km]",colorbar=true,
                                    title="ρ",
                                    aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),                             
                                    ylims=(M.ymin/1e3, M.ymax/1e3),
                                    layout=(3,1),subplot=1)
                        quiver!(p,x.c2d[1:Pl.qinc:end,1:Pl.qinc:end]./1e3,
                            y.c2d[1:Pl.qinc:end,1:Pl.qinc:end]./1e3,
                            quiver=(D.vxc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc,
                                    D.vyc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc),        
                            la=0.5,color="white",layout=(3,1),subplot=1)
                        scatter!(p,Ma.x[1:end]./1e3,Ma.y[1:end]./1e3,
                            ms=3,ma=0.5,mc=Ma.phase[1:end],markerstrokewidth=0.0,
                            xlabel="x[km]",ylabel="y[km]",colorbar=true,
                            title="tracers",label="",
                            xlims=(xwave/1e3-Δ.x/1e3*1.1, xwave/1e3+Δ.x/1e3*1.1), 
                            ylims=(-ywave/1e3-Δ.y/1e3*1.1, -ywave/1e3+Δ.y/1e3*1.1),
                            layout=(3,1),subplot=3)
                        heatmap!(p,x.c./1e3,y.c./1e3,log10.(abs.(D.ηc)'),
                                    color=reverse(cgrad(:roma)),
                                    xlabel="x[km]",ylabel="y[km]",title="η_c",
                                    clims=(15,27),
                                    aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),                             
                                    ylims=(M.ymin/1e3, M.ymax/1e3),colorbar=true,
                                    layout=(2,1),subplot=2)
                        scatter!(p,(xwave/1.0e3,-ywave/1.0e3),
                                markersize=3,label="",color=:black,
                                layout=(3,1),subplot=3)
                        scatter!(p,(x1.vy2d[xn+1,yn]/1e3,y1.vy2d[xn+1,yn]/1e3),
                                markersize=3,label="",color=:blue,
                                layout=(3,1),subplot=3)
                        scatter!(p,(x1.vy2d[xn+2,yn]/1e3,y1.vy2d[xn+2,yn]/1e3),
                                markersize=3,label="",color=:blue,
                                layout=(3,1),subplot=3)
                        scatter!(p,(x1.vy2d[xn+1,yn+1]/1e3,y1.vy2d[xn+1,yn+1]/1e3),
                                markersize=3,label="",color=:red,
                                layout=(3,1),subplot=3)
                        scatter!(p,(x1.vy2d[xn+2,yn+1]/1e3,y1.vy2d[xn+2,yn+1]/1e3),
                                markersize=3,label="",color=:red,
                                layout=(3,1),subplot=3)
                        display(p)
                    end
                    if k == 1 && n == 1 && o == 1
                        scatter!(q,(1/nm[k]/nm[k],PP.ε[1]),
                                    ms=ms[1],markershape=:circle,
                                    label=string(-δA," [m]"),color=mc[l],
                                    xlabel="1/nmx/nmy",ylabel="ε [ % ]",
                                    xscale=:log10, yscale=:log10,
                                    title=string("ηᵣ = ",ηᵣ[o]),
                                    xlims=(1/(maximum(nm)+3)/(maximum(nm)+3), .5),
                                    ylims=(1e-1, 1e2),
                                    layout=(size(addnoise,2),size(ηᵣ,2)),
                                    subplot=((n-1)*size(ηᵣ,2)+o))
                    else
                        scatter!(q,(1/nm[k]/nm[k],PP.ε[1]),
                                    ms=ms[1],markershape=:circle,
                                    color=mc[l],label="",
                                    xlabel="1/nmx/nmy",ylabel="ε [ % ]",
                                    xscale=:log10,yscale=:log10,
                                    title=string("ηᵣ = ",ηᵣ[o]),
                                    xlims=(1/(maximum(nm)+3)/(maximum(nm)+3), .5),
                                    ylims=(1e-1, 1e2),
                                    layout=(size(addnoise,2),size(ηᵣ,2)),
                                    subplot=((n-1)*size(ηᵣ,2)+o))
                    end
                end # Loop δA - l
            end # Loop nm - k 
        end # Loop ηᵣ - o
    end # Loop addnoise - n
    if save_fig == 1
        savefig(q,string("./examples/StokesEquation/2D/Results/RTI_Growth_Rate_Res_Test_const_NC.png"))
    else
        display(q)
    end
end # function

RTI_GrowthRate()
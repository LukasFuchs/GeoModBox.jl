using Plots
# using ExtendableSparse
using GeoModBox
using GeoModBox.InitialCondition #, GeoModBox.MomentumEquation.TwoD
# using GeoModBox.AdvectionEquation.TwoD
using GeoModBox.Tracers.TwoD
using Base.Threads
using Printf, LinearAlgebra

function RTI_TracerField()
    plot_fields     =:no
    # avgm            =[:arith]     # Averaging Method for η - default arith
    # s = plot(layout=(3,2))
    q = plot(layout=(3,2))
    ηᵣ          =   1e-6 #[1e-6 1 500]            #   Viscosity ratio
    for w = 1:3
    if w == 1
        avgm =:arith 
    elseif w==2
        avgm=:geom 
    elseif w==3
        avgm=:harm
    end
    @show avgm
    # Define Initial Condition ========================================== #
    addnoise    =   0 # [0 1]
    # Density Averaging ---
    #   centroids or vertices
    ρavg        =   :centroids  
    nc          =   60 # [10 20 40 60 80 100 120 140 160 180 200]
    # Initial Marker distribution ---
    Ini         =   (p=:RTI,) 
    # Perturbation wavelength [ m ]
    λ           =   4e3
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
    
    phase       =   [0,1]
    # ------------------------------------------------------------------- #
    # Divisional factor of the amplitude following Gerya (2009) --------- #
    delfac      =   15 # [15 150 1500]
    ms          =   zeros(3)
    ms          =   [3]
    # mc          =   ["black","red","yellow"] # δA
    # ------------------------------------------------------------------- #
    # Geometry ========================================================== #
    M       =   Geometry(
        ymin    =   -3.0e3,     #   [ m ]
        ymax    =   0.0,
        xmin    =   0.0,
    )
    ar          =   Int64(round(2 * λ / (M.ymax-M.ymin)))  # aspect ratio
    M.xmax      =   (M.ymax-M.ymin)*ar
    @printf("   xmax: %g \n",M.xmax)
    # ------------------------------------------------------------------- #
    # Grid ====================================================== # 
    NC  =   (
        x   =   nc[1],
        y   =   nc[1],
    )
    NV  =   (
        x   =   NC.x + 1,
        y   =   NC.y + 1,
    )
    Δ       =   GridSpacing(
        x   =   (M.xmax - M.xmin)/NC.x,
        y   =   (M.ymax - M.ymin)/NC.y,
    )
    @printf("    Δx = %g, Δy = %g\n",Δ.x,Δ.y)
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
    # ----------------------------------------------------------- #
    # Allocation ================================================ #
    D       =   (
        ρ       =   zeros(Float64,(NC...)),
        ρc2     =   zeros(Float64,(NC...)),
        ρv      =   zeros(Float64,NV...),
        ρe      =   zeros(Float64,(NC.x+2,NC.y+2)),
        p       =   zeros(Float64,(NC...)),
        wte     =   zeros(Float64,(NC.x+2,NC.y+2)),
        wtv     =   zeros(Float64,(NV.x,NV.y)),
        ηc      =   zeros(Float64,NC...),
        ηce     =   zeros(Float64,(NC.x+2,NC.y+2)),
        ηv      =   zeros(Float64,NV...),
    )
    # ----------------------------------------------------------- #
    # Marker Allocation ========================================= #
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
    # Physics =================================================== #
    @printf("    ηᵣ = %g\n",ηᵣ[1])
    # 0 - upper layer; 1 - lower layer
    η₀      =   η₁*ηᵣ[1]    #   Viscosity composition 0 [ Pa s ]
    η       =   [η₀,η₁]     #   Viscosity for phases 
    @printf("    η₀ = %g\n",η₀)
    # ----------------------------------------------------------- #
    δA          =   -(M.ymax-M.ymin)/2/delfac[1]    #   Amplitude [ m ]
    @printf("    δA = %g\n",δA)
    # Tracer Advection ================================== #
    nmx,nmy =   5,5
    noise   =   addnoise[1]
    nmark   =   nmx*nmy*NC.x*NC.y
    Aparam  =   :phase
    # --------------------------------------------------- #
    # Initialize Tracer Position ---
    Ma      =   IniTracer2D(Aparam,nmx,nmy,Δ,M,NC,noise,Ini.p,phase;λ,δA)
    # Count tracer per cell ---
    CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,NV,1)
    # Interpolate density --- 
    # if ρavg==:centroids
        # Interpolate density from markers to cell ---
        Markers2Cells(Ma,nmark,MAVG.PC_th,D.ρe,MAVG.wte_th,D.wte,x,y,Δ,Aparam,ρ;avgm)
        D.ρ     .=   D.ρe[2:end-1,2:end-1]  
    # elseif ρavg==:vertices 
        # Interpolate density from markers to vertices ---
        Markers2Vertices(Ma,nmark,MAVG.PV_th,D.ρv,MAVG.wtv_th,D.wtv,x,y,Δ,Aparam,ρ;avgm)
        # for i = 1:NC.x
        #     D.ρc2[i,:]    .=   (D.ρv[i,1:end-1] .+ 
        #                         D.ρv[i,2:end,:] .+ 
        #                         D.ρv[i+1,1:end-1] .+ 
        #                         D.ρv[i+1,2:end])./4
        # end
    # end     
    # --------------------------------------------------- #
    # Calculate diapir growth rate ---
    xwave       =   (M.xmax-M.xmin)/2  
    ywave       =   (M.ymax-M.ymin)/2 + δA

    xn          =   Int64(floor((xwave+Δ.x/2)/Δ.x))
    yn          =   Int64(floor(((M.ymax-M.ymin)-ywave)/Δ.y)) + 1

    # dx          =   (xwave+Δ.x/2)/Δ.x - xn
    # dy          =   abs(((M.ymax-M.ymin)-ywave)/Δ.y - yn)

    # wvy     =   (1.0-dx)*(1.0-dy) * D.vy[xn+1,yn] + 
    #                 dx*(1.0-dy) * D.vy[xn+2,yn] + 
    #                 (1.0-dx)*dy * D.vy[xn+1,yn+1] + 
    #                 dx*dy * D.vy[xn+2,yn+1]

    # PP.Q[1] =   (ρ₀-ρ₁)*(M.ymax-M.ymin)/2.0*g/2.0/η₁
    # PP.K[1] =   abs(wvy)/abs(δA)/PP.Q[1]
    # PP.ϕ[1] =   2*π*(M.ymax-M.ymin)/2/λ
    # PP.ε[1] =   abs((PP.Kₐ[o]-PP.K[1])/PP.Kₐ[o])*100.0
    xmin    =   19*(M.xmax-M.xmin)/1e3/40
    xmax    =   21*(M.xmax-M.xmin)/1e3/40
    ymin    =   10*(M.ymin-M.ymax)/1e3/20
    ymax    =   9*(M.ymin-M.ymax)/1e3/20
    heatmap!(q,x.c/1e3,y.c/1e3,D.ρ',layout=(3,2),subplot=(w-1)*2+1)
    scatter!(q,(x.c2d/1e3,y.c2d/1e3),markersize=2,markershape=:circle,
                    markercolor=:black,label="",layout=(3,2),subplot=(w-1)*2+1)
    scatter!(q,(x.v2d/1e3,y.v2d/1e3),markersize=2,markershape=:cross,
                    markercolor=:red,label="",layout=(3,2),subplot=(w-1)*2+1)
    scatter!(q,Ma.x[Ma.phase.==1]./1e3,Ma.y[Ma.phase.==1]./1e3,
            ms=3,ma=0.5,mc=Ma.phase[Ma.phase.==1],markerstrokewidth=0.0,
            colorbar=false,label="",
            xlims=(xmin, xmax),
            ylims=(ymin, ymax),layout=(3,2),subplot=(w-1)*2+1)
    contour!(q,x.c/1e3,y.c/1e3,D.ρ',levels=1,color=:white,
                    xlabel="x[km]",ylabel="y[km]",colorbar=false,
                    title="ρ",
                    aspect_ratio=:equal,
                    xlims=(xmin, xmax),
                    ylims=(ymin, ymax),layout=(3,2),subplot=(w-1)*2+1)
    contour!(q,x.v/1e3,y.v/1e3,D.ρv',levels=1,linestyle=:dash,color=:red, 
                    xlabel="x[km]",ylabel="y[km]",colorbar=false,
                    title=string("ρ, ",avgm),
                    aspect_ratio=:equal,
                    xlims=(xmin, xmax),
                    ylims=(ymin, ymax),layout=(3,2),subplot=(w-1)*2+1)
    # display(q)
    # Interpolate Viscosity ---
    Markers2Cells(Ma,nmark,MAVG.PC_th,D.ηce,MAVG.wte_th,D.wte,x,y,Δ,Aparam,η;avgm)
    D.ηc    .=   D.ηce[2:end-1,2:end-1]
    Markers2Vertices(Ma,nmark,MAVG.PV_th,D.ηv,MAVG.wtv_th,D.wtv,x,y,Δ,Aparam,η;avgm)
    heatmap!(q,x.v/1e3,y.v/1e3,(D.ηv'),layout=(3,2),subplot=(w-1)*2+2)
    scatter!(q,(x.c2d/1e3,y.c2d/1e3),markersize=2,markershape=:circle,
                    markercolor=:black,label="",layout=(3,2),subplot=(w-1)*2+2)
    scatter!(q,(x.v2d/1e3,y.v2d/1e3),markersize=2,markershape=:cross,
                    markercolor=:red,label="",layout=(3,2),subplot=(w-1)*2+2)
    
    scatter!(q,Ma.x[Ma.phase.==1]./1e3,Ma.y[Ma.phase.==1]./1e3,
            ms=3,ma=0.5,mc=Ma.phase[Ma.phase.==1],markerstrokewidth=0.0,
            colorbar=false,label="",
            xlims=(xmin, xmax),
            ylims=(ymin, ymax),layout=(3,2),subplot=(w-1)*2+2)
    contour!(q,x.c/1e3,y.c/1e3,(D.ηc'),levels=1,color=:white, 
            xlabel="x[km]",ylabel="y[km]",colorbar=false,
            title=string("η, ",avgm),
            aspect_ratio=:equal,
            xlims=(xmin, xmax),
            ylims=(ymin, ymax),layout=(3,2),subplot=(w-1)*2+2)
    contour!(q,x.v/1e3,y.v/1e3,(D.ηv'),levels=1,linestyle=:dash,color=:red, 
            # xlabel="x[km]",ylabel="y[km]",colorbar=false,
            # title=string("η, ",avgm),
            aspect_ratio=:equal,
            xlims=(xmin, xmax),
            ylims=(ymin, ymax),layout=(3,2),subplot=(w-1)*2+2)
    scatter!(q,(xwave/1.0e3,-ywave/1.0e3),
            markersize=3,label="",color=:black,layout=(3,2),subplot=(w-1)*2+2)
    scatter!(q,(x1.vy2d[xn+1,yn]/1e3,y1.vy2d[xn+1,yn]/1e3),
            markersize=3,label="",color=:blue,layout=(3,2),subplot=(w-1)*2+2)
    scatter!(q,(x1.vy2d[xn+2,yn]/1e3,y1.vy2d[xn+2,yn]/1e3),
            markersize=3,label="",color=:blue,layout=(3,2),subplot=(w-1)*2+2)
    scatter!(q,(x1.vy2d[xn+1,yn+1]/1e3,y1.vy2d[xn+1,yn+1]/1e3),
            markersize=3,label="",color=:red,layout=(3,2),subplot=(w-1)*2+2)
    scatter!(q,(x1.vy2d[xn+2,yn+1]/1e3,y1.vy2d[xn+2,yn+1]/1e3),
            markersize=3,label="",color=:red,layout=(3,2),subplot=(w-1)*2+2)
    # display(s)
    
    if plot_fields==:yes
        p = heatmap(x.c./1e3,y.c./1e3,D.ρ',color=:inferno,
                    xlabel="x[km]",ylabel="y[km]",colorbar=true,
                    title="ρ",
                    aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),                             
                    ylims=(M.ymin/1e3, M.ymax/1e3),
                    layout=(3,1),subplot=1)
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
    
    end
    @printf("    ηᵣ = %g\n",ηᵣ[1])
    display(q)
    savefig(q,string("./examples/StokesEquation/2D/Results/AverageSchemes",ηᵣ[1],".png"))
end # function

RTI_TracerField()
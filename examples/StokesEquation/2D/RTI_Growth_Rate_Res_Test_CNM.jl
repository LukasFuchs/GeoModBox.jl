using Plots
using ExtendableSparse
using GeoModBox
using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
using GeoModBox.AdvectionEquation.TwoD
using GeoModBox.Tracers.TwoD
using Base.Threads
using Printf, LinearAlgebra
using TimerOutputs
using Random

function UpdateRheo(D,ρ,η,avg)

    # Density --- artihmetic averaging ---
    @. D.ρ     =   ρ[1]*(1.0 - D.p) + ρ[2]*D.p
    D.ρe[2:end-1,2:end-1]     .=  D.ρ
    D.ρe[1,:]     .=  D.ρe[2,:]
    D.ρe[end,:]   .=  D.ρe[end-1,:]
    D.ρe[:,1]     .=  D.ρe[:,2]
    D.ρe[:,end]   .=  D.ρe[:,end-1]

    # Viscosity - Centroids and vertices ---
    if avg ==     :arith
        @. D.ηc =   (1.0 - D.p) * η[1] + D.p * η[2]
    elseif avg == :harm
        @. D.ηc =   1.0 / ( (1.0 - D.p) / η[1] + D.p / η[2] )
    elseif avg == :geom
        @. D.ηc =   η[1]^(1.0 - D.p) * η[2]^D.p
    else
        error("Unknown viscosity averaging: $(avg)")
    end
    # --- Extended Centroids-
    D.ηce[2:end-1,2:end-1]     .=  D.ηc
    D.ηce[1,:]     .=  D.ηce[2,:]
    D.ηce[end,:]   .=  D.ηce[end-1,:]
    D.ηce[:,1]     .=  D.ηce[:,2]
    D.ηce[:,end]   .=  D.ηce[:,end-1]
    # --- Vertices -
    if avg == :arith
        @. D.ηv =
            0.25*(
                D.ηce[1:end-1,1:end-1] + 
                D.ηce[2:end  ,1:end-1] + 
                D.ηce[1:end-1,2:end  ] + 
                D.ηce[2:end  ,2:end  ]
            )
    elseif avg == :harm
        @. D.ηv =
            4.0/(
                1/D.ηce[1:end-1,1:end-1] + 
                1/D.ηce[2:end  ,1:end-1] + 
                1/D.ηce[1:end-1,2:end  ] + 
                1/D.ηce[2:end  ,2:end  ]
            )
    elseif avg == :geom
        @. D.ηv =
            exp(0.25*(
                log(D.ηce[1:end-1,1:end-1]) + 
                log(D.ηce[2:end  ,1:end-1]) + 
                log(D.ηce[1:end-1,2:end  ]) + 
                log(D.ηce[2:end  ,2:end  ])
            ))
    else
        error("Unknown viscosity averaging: $(avg)")
    end

    return D
end

function RTI_GrowthRate()
    to              = TimerOutput()
    @timeit to "Ini" begin
    plot_fields     =:no
    save_fig        = 1
    Pl  =   (
        qinc    =   5, 
        qsc     =   100*(60*60*24*365.25)*5e1,
    )
    # How to update the numerical nodes from the marker ---
    #   1) Phase ratio (:PhaseRatio), or
    #   2) Direct bilinear interpolation of marker properties to the grid
    #      (:MarkerInterpolation)
    MaterialInterpolation   =   :MarkerInterpolation
    # Define the averaging scheme used for both material-transfer methods:
    #   1) Arithmetic   (:arith)
    #   2) Geometric    (:geom)
    #   3) Harmonic     (:harm)
    MaterialInterpolation ∈ (:PhaseRatio, :MarkerInterpolation) ||
        error("Unknown material interpolation method: $MaterialInterpolation")
    # Viscosity averaging schemes included in the summary comparison
    avg_list = (:arith, :geom, :harm)
    # Although both strategies employ bilinear interpolation from the 
    # markers to the numerical grid, they differ in the quantity that is 
    # interpolated. The phase-ratio approach first interpolates the 
    # material phase and subsequently reconstructs the material 
    # properties using a mixing law, whereas the direct interpolation 
    # approach transfers the material properties themselves. Since 
    # viscosity is a nonlinear material property, the two approaches 
    # generally produce different numerical solutions.
    # ------------------------------------------------------------------- #
    # Define Initial Condition ========================================== #
    addnoise    =   [0 1]
    nc          =   [10 20 40 60 80 100 120 140 160 180 200]
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
    ηᵣ          =   [1e-6 1 500]            #   Viscosity ratio
    phase       =   [0,1]
    # Script notation:
    #   phase 0: upper layer  (corresponds to η₁, ρ₁, h₁ in Gerya)
    #   phase 1: lower layer  (corresponds to η₂, ρ₂, h₂ in Gerya)
    # ------------------------------------------------------------------- #
    # Divisional factor of the amplitude following Gerya (2019) --------- #
    delfac      =   [15 150 1500]
    ms          =   zeros(3)
    ms          =   [3]
    mc          =   ["black","red","yellow"] # δA
    # Geometry ========================================================== #
    M       =   Geometry(
        ymin    =   -3.0e3,     #   [ m ]
        ymax    =   0.0,
        xmin    =   0.0,
    )
    # Horizontal model width following Gerya (2019): L = 2λ
    M.xmax      =   M.xmin + 2.0 * λ
    @printf("   xmax: %g \n",M.xmax)
    @printf("   λ = %g\n",λ)
    # ------------------------------------------------------------------- #
    # Analytical Solution =============================================== #
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
    ϕ₁          =   (2*π*((M.ymax-M.ymin)/2))/λ
    ϕ₂          =   ϕ₁
    # ------------------------------------------------------------------- #
    # Summary errors:
    # Dimensions:
    # averaging scheme × marker distribution × viscosity ratio × grid resolution
    ε_summary = fill(
        NaN,
        length(avg_list),
        length(addnoise),
        length(ηᵣ),
        length(nc),
    )
    h_grid = zeros(length(nc))
    end
    @timeit to "Averaging Loop" begin
    for ia in eachindex(avg_list)
        avg = avg_list[ia]
        @printf("\nAveraging scheme: %s\n", string(avg))

        # Detailed on-the-fly plot for the current averaging scheme
        q = plot(
            layout = (length(addnoise), length(ηᵣ)),
            size   = (950, 600),
        )

        @timeit to "Noise Loop" begin
        for n in eachindex(addnoise)
            @timeit to "Resolution Loop" begin
            for k in eachindex(nc)
                @printf("    nc = %g\n",nc[k])
                # Grid ====================================================== # 
                NC  =   (
                    x   =   nc[k],
                    y   =   nc[k],
                )
                NV  =   (
                    x   =   NC.x + 1,
                    y   =   NC.y + 1,
                )
                Δ       =   GridSpacing(
                    x   =   (M.xmax - M.xmin)/NC.x,
                    y   =   (M.ymax - M.ymin)/NC.y,
                )
                h_grid[k] = max(Δ.x, Δ.y)
                @printf("    h = %g\n",h_grid[k])
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
                    ρe      =   zeros(Float64,(NC.x+2,NC.y+2)),
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
                    ηce     =   zeros(Float64,(NC.x+2,NC.y+2)),
                    ηv      =   zeros(Float64,NV...),
            p       =   zeros(Float64,NC...),
            p_ex    =   zeros(Float64,(NC.x+2,NC.y+2)),
                )
                # ----------------------------------------------------------- #
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
                # ----------------------------------------------------------- #
                # Boundary Conditions ======================================= #
                VBC     =   (
                    type    =   (E=:freeslip,W=:freeslip,S=:noslip,N=:noslip),
            val     =   (E=zeros(NV.y),W=zeros(NV.y),S=zeros(NV.x),N=zeros(NV.x),
                            vyS=0.0,vyN=0.0,vxW=0.0,vxE=0.0),
                )
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
                # System of Equations ======================================= #
                # Iterations
                niter       =   50
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
                @timeit to "ηr Loop" begin    
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
                    # ------------------------------------------------------- #
                    @timeit to "δA Loop" begin
                    for l in eachindex(delfac) # Loop over perturbation amplitude
                        δA          =   -(M.ymax-M.ymin)/2/delfac[l]    #   Amplitude [ m ]
                        @printf("    δA = %g\n",δA)
                        # Tracer Advection ================================== #
                        @timeit to "Tracer Ini" begin
                        nmx,nmy =   5,5
                        noise   =   addnoise[n]
                        nmark   =   nmx*nmy*NC.x*NC.y
                        Aparam  =   :phase
                        # --------------------------------------------------- #
                        # Initialize Tracer Position ---
                        if noise == 1
                            # Identical realization for all averaging schemes and viscosity ratios
                            Random.seed!(10_000*k + 100*l + n)
                        end
                        Ma      =   IniTracer2D(Aparam,nmx,nmy,Δ,M,NC,noise,Ini.p,phase;λ,δA)
                        # Count tracer per cell ---
                        CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,NV)
                        if MaterialInterpolation == :PhaseRatio
                            # Phase ---
                            Markers2Cells(Ma,nmark,MAVG.PC_th,D.p_ex,MAVG.wte_th,D.wte,x,y,Δ,Aparam,phase)
                            D.p     .=  D.p_ex[2:end-1,2:end-1]
                            UpdateRheo(D,ρ,η,avg)
                        elseif MaterialInterpolation == :MarkerInterpolation
                            # Interpolate from markers to cell ---
                            Markers2Cells(Ma,nmark,MAVG.PC_th,D.p_ex,MAVG.wte_th,D.wte,x,y,Δ,Aparam,phase)
                            D.p     .=  D.p_ex[2:end-1,2:end-1]
                            Markers2Cells(Ma,nmark,MAVG.PC_th,D.ρe,MAVG.wte_th,D.wte,x,y,Δ,Aparam,ρ)
                            D.ρ     .=   D.ρe[2:end-1,2:end-1]  
                            Markers2Cells(Ma,nmark,MAVG.PC_th,D.ηce,MAVG.wte_th,D.wte,x,y,Δ,Aparam,η;avgm=avg)
                        D.ηc    .=   D.ηce[2:end-1,2:end-1]
                            Markers2Vertices(Ma,nmark,MAVG.PV_th,D.ηv,MAVG.wtv_th,D.wtv,x,y,Δ,Aparam,η;avgm=avg)
                        else
                            error("Unknown material interpolation method: $MaterialInterpolation")
                        end
                        end
                        # --------------------------------------------------- #
                        # --------------------------------------------------- #
                        # Momentum Equation ===
                        D.vx    .=  0.0
                        D.vy    .=  0.0
                        D.Pt    .=  1.0
                        @timeit to "Solution Iteration" begin
                        # Assemble Coefficients ========================================= #
                        @timeit to "Assembly" begin
                        K       =   Assembly(NC, NV, Δ, D.ηc, D.ηv, VBC, Num)
                        Kfac    =   lu(K.cscmatrix)
                        end
                        for iter=1:niter
                            # Initial Residual ------------------------------ #
                            @timeit to "Residual" begin
                            Residuals2D!(D,VBC,ε,τ,divV,Δ,D.ηc,D.ηv,g,Fm,FPt)
                            F[Num.Vx]   .=   Fm.x
                            F[Num.Vy]   .=   Fm.y
                            F[Num.Pt]   .=   FPt
                            @printf("||R|| = %1.4e\n", norm(F)/length(F))
                            norm(F)/length(F) < ϵ ? break : nothing
                            end
                            # Solution of the linear system ========================= #
                            @timeit to "Solution" begin
                            δx      =   - (Kfac \ F)
                            end
                            # ----------------------------------------------- #
                            # Update Unknown Variables ====================== #
                            D.vx[:,2:end-1]     .+=  δx[Num.Vx]
                            D.vy[2:end-1,:]     .+=  δx[Num.Vy]
                            D.Pt                .+=  δx[Num.Pt]
                        end
                        end
                        # --------------------------------------------------- #
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
                        @timeit to "Calculate GR" begin
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
                        # Store only:
                        # Store the smallest perturbation for both marker distributions
                        if l == length(delfac)
                            ε_summary[ia, n, o, k] = PP.ε[1]
                        end
                        end
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
                            scatter!(q,(h_grid[k],PP.ε[1]),
                                        ms=ms[1],markershape=:circle,
                                        label=string(-δA," [m]"),color=mc[l],
                                        ylabel="ε [ % ]",
                                        xscale=:log10, yscale=:log10,
                                        xlabel = "Grid spacing h [m]",
                                        xlims  = (30.0, 1000.0),
                                        ylims = (1e-2, 1e7),
                                        layout=(size(addnoise,2),size(ηᵣ,2)),
                                        subplot=((n-1)*size(ηᵣ,2)+o))
                        else
                            scatter!(q,(h_grid[k],PP.ε[1]),
                                        ms=ms[1],markershape=:circle,
                                        color=mc[l],label="",
                                        ylabel="ε [ % ]",
                                        xscale=:log10,yscale=:log10,
                                        xlabel = "Grid spacing h [m]",
                                        xlims  = (30.0, 1000.0),
                                        ylims = (1e-2, 1e7),
                                        layout=(size(addnoise,2),size(ηᵣ,2)),
                                        subplot=((n-1)*size(ηᵣ,2)+o))
                        end
                    end # Loop δA - l
                    end
                end # Loop ηᵣ - o
                end
            end # Loop nc - k
            end
        end # Loop addnoise - n
        end
        # Row labels -------------------------------------------------------- #
        for sp in 1:(length(addnoise) * length(ηᵣ))
            hline!(q[sp], [1.0],  linestyle=:dash,    label="")
            hline!(q[sp], [5.0],  linestyle=:dot,     label="")
            hline!(q[sp], [10.0], linestyle=:dashdot, label="")
        end
        # Column titles
        for o in eachindex(ηᵣ)
            plot!(q[o], title="ηᵣ = $(ηᵣ[o])")
        end
        plot!(q[1], ylabel="Regular markers\n\nε [%]")
        plot!(q[4], ylabel="Perturbed markers\n\nε [%]")

        for sp in (2, 3, 5, 6)
            plot!(q[sp], ylabel="")
        end
        for sp in (1, 2, 3)
            plot!(q[sp], xlabel="")
        end
        # ------------------------------------------------------------------ #
    if save_fig == 1
        savefig(q,string("./examples/StokesEquation/2D/Results/RTI_Growth_Rate_Res_Test_const_NM_",
                    avg,"_",MaterialInterpolation,".png"))
    else
        display(q)
    end
        end # Averaging loop
    end # Averaging timing block
    # Summary plot ====================================================== #
    grid_spacing = copy(h_grid)

    avg_shapes = (
        arith = :circle,
        geom  = :diamond,
        harm  = :utriangle,
    )

    avg_labels = (
        arith = "Arithmetic",
        geom  = "Geometric",
        harm  = "Harmonic",
    )

    q_summary = plot(
        layout        = (length(addnoise), length(ηᵣ)),
        size          = (1100, 700),
        left_margin   = 8Plots.mm,
        bottom_margin = 7Plots.mm,
        top_margin    = 4Plots.mm,
        legendfontsize = 8,
    )

    for n in eachindex(addnoise)
        for o in eachindex(ηᵣ)

            subplot_index = (n - 1) * length(ηᵣ) + o

            for ia in eachindex(avg_list)

                avg_current = avg_list[ia]

                # Force x and y to be one-dimensional vectors.
                # Each averaging scheme is therefore plotted as one series.
                xdata = vec(grid_spacing)
                ydata = vec(ε_summary[ia, n, o, :])

                plot!(
                    q_summary[subplot_index],
                    xdata,
                    ydata,
                    marker      = avg_shapes[avg_current],
                    markersize  = 6,
                    linewidth   = 1.5,
                    linestyle   = :solid,
                    xscale      = :log10,
                    yscale      = :log10,
                    xlims       = (30.0, 1000.0),
                    ylims       = (1e-2, 1e7),

                    # Show the legend only once, in the upper-left panel.
                    label        = n == 1 && o == 1 ?
                                avg_labels[avg_current] : "",
                    legend       = n == 1 && o == 1 ?
                                :topright : false,

                    xlabel       = n == length(addnoise) ?
                                "Grid spacing h [m]" : "",
                    ylabel       = o == 1 ?
                                "Relative error ε [%]" : "",
                    title        = n == 1 ?
                                "ηᵣ = $(ηᵣ[o])" : "",
                )
            end

            hline!(
                q_summary[subplot_index],
                [1.0],
                linestyle = :dash,
                linewidth = 1,
                label     = "",
            )

            hline!(
                q_summary[subplot_index],
                [10.0],
                linestyle = :dot,
                linewidth = 1,
                label     = "",
            )
        end
    end

    # Row labels
    plot!(
        q_summary[1],
        ylabel = "Regular markers\n\nRelative error ε [%]",
    )

    plot!(
        q_summary[length(ηᵣ) + 1],
        ylabel = "Perturbed markers\n\nRelative error ε [%]",
    )

    # Use identical limits in all six panels
    for sp in 1:(length(addnoise) * length(ηᵣ))
        plot!(
            q_summary[sp],
            xscale = :log10,
            yscale = :log10,
            xlims = (30.0, 1000.0),
            ylims  = (1e-2, 1e7),
        )
    end

    if save_fig == 1
        savefig(
            q_summary,
            string(
                "./examples/StokesEquation/2D/Results/",
                "RTI_Growth_Rate_Averaging_Summary_const_NM_",
                MaterialInterpolation,
                ".png",
            ),
        )
    else
        display(q_summary)
    end
    # ------------------------------------------------------------------ #
    display(to)
end # function

RTI_GrowthRate()
using Plots, ExtendableSparse, Base.Threads
using GeoModBox, GeoModBox.Tracers.TwoD
using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD
using GeoModBox.AdvectionEquation.TwoD, GeoModBox.HeatEquation.TwoD
using Statistics
using Printf, LinearAlgebra, TimerOutputs, Interpolations, LsqFit

# ======================================================================= #
# ======================= HELPER FUNCTIONS ============================== #
# ======================================================================= #
# ----------------------------------------------------------------------- #
function GetSecondInvariants!(ε,τ)
    # Get shear strain rate and stress on centroids ---
    @. ε.xyc    =   (ε.xy[1:end-1,1:end-1] + ε.xy[2:end-0,1:end-1] + 
                    ε.xy[1:end-1,2:end-0] + ε.xy[2:end-0,2:end-0])/4.0
    @. τ.xyc    =   (τ.xy[1:end-1,1:end-1] + τ.xy[2:end-0,1:end-1] + 
                    τ.xy[1:end-1,2:end-0] + τ.xy[2:end-0,2:end-0])/4.0
    # Get invariants ---
    @. ε.II     =   sqrt(0.5*(ε.xx^2 + ε.yy^2) + ε.xyc^2)
    @. τ.II     =   sqrt(0.5*(τ.xx^2 + τ.yy^2) + τ.xyc^2)
end
# ----------------------------------------------------------------------- #
function GetTimeStep!(T,Δ,P,D)
    T.Δc    =   T.Δfacc * minimum((Δ.x,Δ.y)) / 
                    (sqrt(maximum(abs.(D.vx))^2 + maximum(abs.(D.vy))^2))
    T.Δd    =   T.Δfacd * (1.0 / (2.0 * P.κ *(1.0/Δ.x^2 + 1/Δ.y^2)))
    T.Δ     =   minimum([T.Δd,T.Δc])
end
# ----------------------------------------------------------------------- #
function ComputeStrainStress2D!(D,ε,τ,divV,Δ)
    @. divV =   (D.vx[2:end,2:end-1] - D.vx[1:end-1,2:end-1])/Δ.x + (D.vy[2:end-1,2:end] - D.vy[2:end-1,1:end-1])/Δ.y
    @. ε.xx =   (D.vx[2:end,2:end-1] - D.vx[1:end-1,2:end-1])/Δ.x - 1.0/3.0*divV
    @. ε.yy =   (D.vy[2:end-1,2:end] - D.vy[2:end-1,1:end-1])/Δ.y - 1.0/3.0*divV
    @. ε.xy =   0.5*( (D.vx[:,2:end] - D.vx[:,1:end-1])/Δ.y + (D.vy[2:end,:] - D.vy[1:end-1,:])/Δ.x ) 
    @. τ.xx =   2.0 * D.ηc * ε.xx
    @. τ.yy =   2.0 * D.ηc * ε.yy
    @. τ.xy =   2.0 * D.ηv * ε.xy
end
# ----------------------------------------------------------------------- #
function UpdateRheo!(ε,Rhe,D,P;ini=:no)

    # Set effective strain rate to avoid Infs and NaNs --- 
    @. Rhe.εIIeff   =   max(ε.II, Rhe.εIImin)

    # Update rheology ---
    # Viscosity ---            
    @. Rhe.ηmat =   1e20 #(2.0^((1.0-Rhe.n[1])/Rhe.n[1]))/(3.0^((1.0+Rhe.n[1])/(2.0*Rhe.n[1]))) * 
                    #  Rhe.A[1]^(-1/Rhe.n[1]) * Rhe.εIIeff^(1/Rhe.n[1]-1) * exp(Rhe.Ea[1]/(Rhe.n[1]*P.RG*(D.T+273.15)))
    @. Rhe.ηinc =   1e20 # (2.0^((1.0-Rhe.n[2])/Rhe.n[2]))/(3.0^((1.0+Rhe.n[2])/(2.0*Rhe.n[2]))) * 
                    #  Rhe.A[2]^(-1/Rhe.n[2]) * Rhe.εIIeff^(1/Rhe.n[2]-1) * exp(Rhe.Ea[2]/(Rhe.n[2]*P.RG*(D.T+273.15)))
    
    # @. Rhe.ηnew =   Rhe.ηmat.^(1.0 .- D.p) .* Rhe.ηinc.^D.p
    if Rhe.avg_p == :arithmetic
        @. Rhe.ηnew = (1.0 - D.p) * Rhe.ηmat + D.p * Rhe.ηinc
    elseif Rhe.avg_p == :harmonic
        @. Rhe.ηnew = 1.0 / ( (1.0 - D.p) / Rhe.ηmat + D.p / Rhe.ηinc )
    elseif Rhe.avg_p == :geometric
        @. Rhe.ηnew = Rhe.ηmat^(1.0 - D.p) * Rhe.ηinc^D.p
    else
        error("Unknown viscosity averaging: $(Rhe.avg_p)")
    end

    if ini==:no
        @. D.ηc     =   exp((1-Rhe.ω)*log(D.ηc) + Rhe.ω*log(Rhe.ηnew))
    else
        @. D.ηc     =   Rhe.ηnew
    end
    # --- Extended Centroids-
    D.η_ex[2:end-1,2:end-1]     .=  D.ηc
    D.η_ex[1,:]     .=  D.η_ex[2,:]
    D.η_ex[end,:]   .=  D.η_ex[end-1,:]
    D.η_ex[:,1]     .=  D.η_ex[:,2]
    D.η_ex[:,end]   .=  D.η_ex[:,end-1]
    # --- Vertices -
    if Rhe.avg_v == :arithmetic
        @. D.ηv =
            0.25*(
                D.η_ex[1:end-1,1:end-1] + 
                D.η_ex[2:end  ,1:end-1] + 
                D.η_ex[1:end-1,2:end  ] + 
                D.η_ex[2:end  ,2:end  ]
            )
    elseif Rhe.avg_v == :harmonic
        @. D.ηv =
            4.0/(
                1/D.η_ex[1:end-1,1:end-1] + 
                1/D.η_ex[2:end  ,1:end-1] + 
                1/D.η_ex[1:end-1,2:end  ] + 
                1/D.η_ex[2:end  ,2:end  ]
            )
    elseif Rhe.avg_v == :geometric
        @. D.ηv =
            exp(0.25*(
                log(D.η_ex[1:end-1,1:end-1]) + 
                log(D.η_ex[2:end  ,1:end-1]) + 
                log(D.η_ex[1:end-1,2:end  ]) + 
                log(D.η_ex[2:end  ,2:end  ])
            ))
    else
        error("Unknown viscosity averaging: $(Rhe.avg_v)")
    end
end
# ----------------------------------------------------------------------- #
function Diagnostics!(phbd,halfwidth,nprof,
                        M,T,D,Ini,x,y,ε,
                        θsb,θsb2,Dsb,εf,δTemp,shortening,strain,style,
                        xp,yp,it)
    pp  =   nothing
    r   =   nothing
    # Estimate shear band angle ----------------------------------------- #
    # Matrix shear-zone candidate points
    mask_mat    =   D.p .< phbd
    mask_right  =   x.c2d .> 1.5e3
    mask_top    =   y.c2d .> 1.5e3
    mask_high   =   ε.II .>  Ini.ε
    # Shear band mask --- Total ---
    mask_shear_band = mask_mat .& mask_right .& mask_top .& mask_high
    # Remove points of lower strain rate ---
    mask_low  = ε.II .>  0.4 * maximum(ε.II[mask_shear_band])
    # Final mask ---
    mask_band = mask_shear_band .& mask_low
    # Define profile through the shear band ----------------------------- #
    if style==:moving
        # Shear band point coordinates ---
        xb      =   x.c2d[mask_band]
        yb      =   y.c2d[mask_band]
        # Fit y = a*x + b ---
        Afit    =   hcat(xb, ones(length(xb)))
        a, b    =   Afit \ yb
        # Define center point for profile ------------------------------- #
        λ           =   0.2   # 0 near inclusion, 1 near top shear band points
        xc          =   λ * maximum(xb)
        yc          =   a * xc + b
        θsb[it]     =   atan(a)   # shear band angle
    elseif style==:fixed
        xc          =   3.8e3 * exp(-strain[it])
        yc          =   3.8e3 * exp( strain[it])
        θsb[it]     =   atan(yc, xc)
    end
    # Calculate points for profile ---
    # Normal to the shear band ---
    nx      =   -sin(θsb[it])
    ny      =   cos(θsb[it])
    # Profile length ---
    sbl     =   range(-halfwidth, halfwidth, length=nprof)
    # Profile coordinates ---
    @. xp   =   xc + sbl * nx
    @. yp   =   yc + sbl * ny
    # --------------------------------------------------------------- #
    # Interpolate strain rate and temperature ----------------------- #
    itpε    =   extrapolate(interpolate((x.c, y.c), ε.II, Gridded(Linear())), NaN)
    itpT    =   extrapolate(interpolate((x.c, y.c), D.T,  Gridded(Linear())), NaN)
    εprof   =   [itpε(xp[k], yp[k]) for k in eachindex(xp)]
    Tprof   =   [itpT(xp[k], yp[k]) for k in eachindex(xp)]
    # --------------------------------------------------------------- #
    # Fit gaussian to strain rate profile --------------------------- #
    # Scale strain rate profile --- gives a better fit -
    εscale      =   Ini.ε
    εprof_s     =   εprof ./ εscale

    εbg         =   quantile(εprof_s, 0.1)
    εmax        =   maximum(εprof_s)
    imax        =   argmax(εprof_s)
    s0_guess    =   sbl[imax]
    if style==:fixed
        θsb2[it]    =   atan(yp[imax],xp[imax])
    end

    fitmask = abs.(sbl .- s0_guess) .< 3.0e3

    model(s,p) = p[1] .+ p[2] .* exp.(-((s .- p[3]).^2) ./ (2*p[4]^2))

    p0      =   [εbg, εmax - εbg, s0_guess, 700.0]
    fit     =   curve_fit(model, sbl[fitmask], εprof_s[fitmask], p0)
    pars    =   coef(fit)

    σ       =   abs(pars[4])
    Dsb[it] =   2*σ

    # if style==:fixed
    εf[it]      =   maximum(εprof_s)
    δTemp[it]   =   maximum(Tprof) - minimum(D.T)
    # else
    #     εf[it]      =   maximum(ε.II[mask_band]) ./ εscale
    #     δTemp[it]   =   maximum(D.T[mask_band]) - minimum(D.T)
    # end

    # δTemp[it]   =   maximum(Tprof) - Ini.Tbg
    # δTemp[it]   =   maximum(Tprof) - minimum(D.T)
    # δTemp[it]   =   maximum(D.T) - Ini.Tbg
    # δTemp[it]   =   maximum(D.T[mask_band]) - minimum(D.T)
    # ------------------------------------------------------------------- #
    if mod(it,5) == 0 || it == 1 || it == T.itmax
        # r = heatmap(x.c./1e3,y.c./1e3,log10.(ε.II'),color=cgrad(:batlow),
        #             xlabel="x[km]",ylabel="y[km]",
        #             clims=(-13.5,-12.0),
        #             aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3), 
        #             ylims=(M.ymin/1e3, M.ymax/1e3),colorbar=true,
        #             size = (900,700), titlefontsize = 20,
        #             guidefontsize = 20, tickfontsize = 16,
        #             colorbar_tickfontsize = 14,
        #             colorbar_titlefontsize = 16)
        r = heatmap(x.c./1e3,y.c./1e3,(D.T'),color=cgrad(:lajolla),
                    xlabel="x[km]",ylabel="y[km]",
                    # clims=(-13.5,-12.0),
                    aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3), 
                    ylims=(M.ymin/1e3, M.ymax/1e3),colorbar=true,
                    size = (900,700), titlefontsize = 20,
                    guidefontsize = 20, tickfontsize = 16,
                    colorbar_tickfontsize = 14,
                    colorbar_titlefontsize = 16)
        scatter!(r,[xc/1e3],[yc/1e3],
                    ms=4,ma=0.5,mc=:black,markerstrokewidth=0.0,label="")
        if style==:moving
            scatter!(r,[xb./1e3],[yb./1e3],
                        ms=1,ma=0.5,mc=:yellow,markerstrokewidth=0.0,label="")
        end
        contour!(r,x.c./1e3,y.c./1e3,log10.(ε.II)',levels = [log10(Ini.ε)],
                    color=:white,linewidth=1,la=0.5)
        contour!(r,x.c./1e3,y.c./1e3,((D.p.+0.5).*log10.(Ini.ε))',
                    levels = [log10(Ini.ε)],
                    color=:black,linewidth=2)
        plot!(r,xp./1e3,yp./1e3,color=:red,label="",linewidth=2)
                    
        pp = plot(sbl./1e3,(εprof),
                    xlabel = "s [km]",ylabel = "ε_II",
                    label = "ε_{II}",title="Strain Rate",
                    xlims=(-halfwidth./1e3,halfwidth./1e3),
                    ylims=(1e-15,1e-12),
                    layout=(1,2),subplot=1)
        plot!(pp,sbl./1e3,Tprof,xlabel = "s [km]",
                    ylabel = "T",title="Temperature",
                    xlims=(-halfwidth./1e3,halfwidth./1e3),
                    ylims=(400.0,700.0),
                    label = "",layout=(1,2),subplot=2)
        plot!(pp,sbl./1e3,εscale .* model(sbl,pars),
                lw=1,la=0.5,color=:red,
                label="Gaussian fit",layout=(1,2),subplot=1)
    end
    return pp, r
end
# ======================================================================= #
# ======================================================================= #

# ======================================================================= #
# ========================== MAIN FUNCTION ============================== #
# ======================================================================= #
function ShearHeatingShearBands(Diff,θ,Adv,style)
    to          =   TimerOutput()
    @timeit to "Ini" begin
    # Define Initial Condition ========================================== #
    Ini         =   (T      = :gaussian,               # Temperature
                     Tbg    = 400.0,                # [ °C ]
                     p      = :ShearBandSetting,    # Phasedistribution
                     radius = 3.0e3,                # [ m ]
                     V      = :ShearBandPS,         # Velocity
                     ε      = 5e-14,                # Background strain rate
    ) 
    # ------------------------------------------------------------------- #
    # Define numerical methods ========================================== #
    # if FD.Method.Diff =: dc; θ-rule
    #       θ = 0   -> implicit
    #       θ = 0.5 -> CN discretization
    #       θ = 1.0 -> explicit
    FD  =   (Method = (
                Diff = Diff,
                Adv  = Adv,
                θ    = θ), 
    )
    # Diagnostics ----
    # style   =:fixed
    # ------------------------------------------------------------------- #
    # Geometry ========================================================== #
    M       =   Geometry(
        xmin    =   -35.0e3,    # [ m ]
        xmax    =   35.0e3,     # [ m ]
        ymin    =   0.0e3,      # [ m ]
        ymax    =   40.0e3,     # [ m ]
    )
    # -------------------------------------------------------------------- #
    # Grid =============================================================== #
    NC      =   (
        x   =   200, 
        y   =   100,
    )
    NV      =   (
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
    # -------------------------------------------------------------------- #
    # Physics ============================================================ #
    P   =   Physics(
        g       =   0.0,            #   Gravitational acceleration [m/s^2]
        ρ₀      =   2700.0,         #   Reference density [kg/m^3]
        k       =   2.5,            #   Thermal conductivity [ W/m/K ]
        cp      =   1050.0,         #   Heat capacity [ J/kg/K ]
    )
    P.κ =  0
    # ------------------------------------------------------------------- #
    # Define rheology paramters ========================================= #
    Rhe     =   ( 
        #           [Matrix Inclusion]
        A       =   [3.2e-20 3.16e-26],     #   Pre-exponential factor
        n       =   [3.0 3.3],              #   Stress exponent
        Ea      =   [276e3 186e3],          #   Activation energy
        # Viscosity damping factor ---
        ω       =   0.0,
        # Lower cut off strain rate --- 
        εIImin  =   1e-20,
        # Initialize some arrays ---
        εIIeff  =   zeros(Float64,NC...),
        ηmat    =   zeros(Float64,NC...),
        ηinc    =   zeros(Float64,NC...),
        ηnew    =   zeros(Float64,NC...),
        avg_p   =   :geometric, 
        avg_v   =   :geometric,
    )
    # ------------------------------------------------------------------- #
    # Define phase ID =================================================== #
    phase       =   [0,1]
    # ------------------------------------------------------------------- #
    # Animation and Plot Settings ======================================= #
    save_fig    =   1

    path        =   string("./examples/ShearHeating/2D/Results/",
                            FD.Method.Diff,"_",FD.Method.θ,"_",
                            FD.Method.Adv,"_", NC.x,"_",
                            NC.y,"_",Rhe.avg_p,"_",Rhe.avg_v,"_",
                            style,"_test") #
    if save_fig == 1
        isdir(path) || mkpath(path)
        framepath2D   = joinpath(path, "frames_2D")
        framepathProf = joinpath(path, "frames_profiles")
        isdir(framepath2D)   || mkpath(framepath2D)
        isdir(framepathProf) || mkpath(framepathProf)
        anim2D      =   Plots.Animation(framepath2D, String[])
        animProf    =   Plots.Animation(framepathProf, String[])
        filename    =   string("/ShearHeatingBands")
        dfilen      =   string("Diagnostics_",FD.Method.Diff,"_",
                            FD.Method.θ,"_",FD.Method.Adv,"_",
                            NC.x,"_",NC.y,"_.txt")
        dfile       =   open("$path/$dfilen","w")
        println(dfile,"# it strain shortening ε_f δT D θ θc")
    end
    # ------------------------------------------------------------------- #
    # Allocation ======================================================== #
    D   =   DataFields(
        T       =   zeros(Float64,(NC...)),
        Hs      =   zeros(Float64,NC...),
        T0      =   zeros(Float64,(NC...)),
        T_exD0  =   zeros(Float64,(NC.x+2,NC.y+2)),
        T_ex    =   zeros(Float64,(NC.x+2,NC.y+2)),
        T_ex0   =   zeros(Float64,(NC.x+2,NC.y+2)),
        vx      =   zeros(Float64,NV.x,NC.y+2),
        vy      =   zeros(Float64,NC.x+2,NV.y),
        Pt      =   zeros(Float64,NC...),
        p       =   zeros(Float64,NC...),
        p_ex    =   zeros(Float64,(NC.x+2,NC.y+2)),
        ρ       =   zeros(Float64,NC...),
        ρ_ex    =   zeros(Float64,(NC.x+2,NC.y+2)),
        vxc     =   zeros(Float64,NC...),
        vyc     =   zeros(Float64,NC...),
        vxco    =   zeros(Float64,(NC...)),
        vyco    =   zeros(Float64,(NC...)),
        vc      =   zeros(Float64,NC...),
        wt      =   zeros(Float64,(NC.x,NC.y)),
        wte     =   zeros(Float64,(NC.x+2,NC.y+2)),
        wtv     =   zeros(Float64,(NV.x,NV.y)),
        ηc      =   zeros(Float64,NC...),
        η_ex    =   zeros(Float64,(NC.x+2,NC.y+2)),
        ηv      =   zeros(Float64,NV...),
    )
    # Needed for Euler advection schemes ---
    vx_bg_c     =   zeros(Float64,NC...)
    vy_bg_c     =   zeros(Float64,NC...)
    vxc_res     =   zeros(Float64,NC...)
    vyc_res     =   zeros(Float64,NC...)
    vxco_res    =   zeros(Float64,NC...)
    vyco_res    =   zeros(Float64,NC...)
    # Needed for the defect correction solution ---
    divV        =   zeros(Float64,NC...)
    ε           =   (
        xx      =   zeros(Float64,NC...), 
        yy      =   zeros(Float64,NC...), 
        xy      =   zeros(Float64,NV...),
        xyc     =   zeros(Float64,NC...),
        II      =   zeros(Float64,NC...),
    )
    τ           =   (
        xx      =   zeros(Float64,NC...), 
        yy      =   zeros(Float64,NC...), 
        xy      =   zeros(Float64,NV...),
        xyc     =   zeros(Float64,NC...),
        II      =   zeros(Float64,NC...),
    )
    # ------------------------------------------------------------------- #
    # Boundary Conditions =============================================== #
    VBC     =   (
        type    =   (E=:ps,W=:ps,S=:freeslip,N=:ps),
        val     =   (E=zeros(NV.y),W=zeros(NV.y),S=zeros(NV.x),N=zeros(NV.x),
                    vxE=zeros(NC.y),vxW=zeros(NC.y),vyS=zeros(NC.x),vyN=zeros(NC.x)),
    )
    TBC     =   (
        type    =   (E=:Neumann,W=:Neumann,S=:Neumann,N=:Neumann),
        val     =   (E=zeros(NC.y),W=zeros(NC.y),S=zeros(NC.x),N=zeros(NC.x)),
    )
    # ------------------------------------------------------------------- #
    # Initial Condition ================================================= #
    # Velocity ---
    IniVelocity!(Ini.V,D,VBC,NV,Δ,M,x,y;Ini.ε)
    # Temperature --- 
    IniTemperature!(Ini.T,M,NC,D,x,y;Tb=Ini.Tbg,Ta=Ini.Tbg+200)
    @. ε.II     =   Ini.ε
    # ------------------------------------------------------------------- #
    # Time ============================================================== #
    T   =   TimeParameter( 
        Δfacc   =   0.9,                #   Courant time factor
        Δfacd   =   0.9,                #   Diffusion time factor
        itmax   =   200,                #   Maximum iterations
    )
    Time    =   zeros(T.itmax)
    # Initialize Time Step ---
    GetTimeStep!(T,Δ,P,D)
    # ------------------------------------------------------------------- #
    # Diagnostics ======================================================= #
    εf          =   zeros(T.itmax)
    δTemp       =   zeros(T.itmax)
    θsb         =   zeros(T.itmax)
    θsb2        =   zeros(T.itmax)
    Dsb         =   zeros(T.itmax)
    strain      =   zeros(T.itmax)
    shortening  =   zeros(T.itmax)
    # Setup parameters -------------------------------------------------- #
    phbd        =   0.01    # Phase boundary value
    halfwidth   =   4e3     # Half of the profile length
    nprof       =   200     # Number of points in the profile
    # Profile coordinates ---
    xp          =   zeros(nprof)
    yp          =   zeros(nprof)
    #-------------------------------------------------------------------- #
    end
    # Tracer Initialization ============================================= #
    @timeit to "Tracer Ini" begin
    nmx,nmy     =   3,3
    noise       =   0
    nmark       =   nmx*nmy*NC.x*NC.y
    Aparam      =   :phase
    if Aparam==:thermal && FD.Method.Adv !=:tracers
        @printf("ERROR! Tracer advection needed to transport temperature!")
        return
    end
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
                        ellA=Ini.radius)
    # RK4 weights ---
    rkw     =   1.0/6.0*[1.0 2.0 2.0 1.0]   # for averaging
    rkv     =   1.0/2.0*[1.0 1.0 2.0 2.0]   # for time stepping
    # Count marker per cell ---
    CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,NV)
    # Temperature --- if FD.Method.Adv==:tracers 
    #   -> FromCtoM
    # ------------------------------------------------------------------- #
    # Update cell centroids and vertices (here only density and phase) -- #
    # Update Centroids --- 
    # Phase ---
    Markers2Cells(Ma,nmark,MAVG.PC_th,D.p_ex,MAVG.wte_th,D.wte,x,y,Δ,Aparam,phase)
    D.p     .=  D.p_ex[2:end-1,2:end-1]
    # Density ---
    @. D.ρ     =   P.ρ₀*(1.0 - D.p) + P.ρ₀*D.p
    D.ρ_ex[2:end-1,2:end-1]     .=  D.ρ
    D.ρ_ex[1,:]     .=  D.ρ_ex[2,:]
    D.ρ_ex[end,:]   .=  D.ρ_ex[end-1,:]
    D.ρ_ex[:,1]     .=  D.ρ_ex[:,2]
    D.ρ_ex[:,end]   .=  D.ρ_ex[:,end-1]
    # Viscosity - Centroids and vertices ---
    UpdateRheo!(ε,Rhe,D,P;ini=:yes)
    # ------------------------------------------------------------------- #
    end
    # System of Equations =============================================== #
    # Numbering, without ghost nodes! ---
    off    = [  NV.x*NC.y,                          # vx
                NV.x*NC.y + NC.x*NV.y,              # vy
                NV.x*NC.y + NC.x*NV.y + NC.x*NC.y]  # Pt
    Num    =    (
        Vx  =   reshape(1:NV.x*NC.y, NV.x, NC.y), 
        Vy  =   reshape(off[1]+1:off[1]+NC.x*NV.y, NC.x, NV.y), 
        Pt  =   reshape(off[2]+1:off[2]+NC.x*NC.y,NC...),#
        T   =   reshape(1:NC.x*NC.y, NC.x, NC.y),
    )
    ndof    =   maximum(Num.T)  
    # Momentum conservation (MCE) ---
    # Iterations --- 
    niterM      =   50          #   Number of iterations 
    atolM       =   1e-8        #   Absolute tolerance
    rtolM       =   1e-5        #   Relative tolerance; r = atolM0/atolM
    RM          =   0.0         #   Initialize absolute residual    
    RMrel       =   0.0         #   Initialize relative residual 
    KM          =   ExtendableSparseMatrix(ndof,ndof)
    # Residuals ---
    Fm     =    (
        x       =   zeros(Float64,NV.x, NC.y), 
        y       =   zeros(Float64,NC.x, NV.y)
    )
    FPt     =   zeros(Float64,NC...)
    # Energy conservation (ECE) ---
    if FD.Method.Diff==:implicit || FD.Method.Diff==:CN
        if FD.Method.Diff==:CN
            K1      =   ExtendableSparseMatrix(ndof,ndof)
            K2      =   ExtendableSparseMatrix(ndof,ndof)
        else
            K       =   ExtendableSparseMatrix(ndof,ndof)
        end
        rhs         =   zeros(ndof)
    elseif FD.Method.Diff==:dc
        niter       =   10
        ϵ           =   1e-20
        K           =   ExtendableSparseMatrix(ndof,ndof)
        R           =   zeros(Float64,NC...)
        ∂2T         =   (∂x2=zeros(NC.x, NC.y), ∂y2=zeros(NC.x, NC.y),
                        ∂x20=zeros(NC.x, NC.y), ∂y20=zeros(NC.x, NC.y))
    end
    # ------------------------------------------------------------------- #
    # Time Loop ========================================================= #
    @timeit to "Time Loop" begin
    for it = 1:T.itmax
        δx      =   zeros(maximum(Num.Pt))
        F       =   zeros(maximum(Num.Pt))
        R0      =   0.0
        # Update Time ---
        if it>1
            Time[it]    =   Time[it-1] + T.Δ    
            strain[it]  =   strain[it-1] + Ini.ε * T.Δ
        end
        shortening[it]  =   100 .* (1 .- exp.(-strain[it]))
        # ---
        @printf("\nTime step: #%04d, Time [Myr]: %04e, Shortening: %2.2f\n ",it,
                    Time[it]/(60*60*24*365.25)/1.0e6,shortening[it])
        # Momentum Conservation Equation ================================ #
        # Initial Residual ---
        if it == 1
            D.vx[2:end-1,:]    .=  0.0
            D.vy[:,1:end-1]    .=  0.0
            D.Pt               .=  0.0
        end
        @timeit to "Solution Iteration" begin
        @printf("---Momentum Calculation ---\n")        
        for iter = 1:niterM
            @timeit to "Residual" begin
            Residuals2D!(D,VBC,ε,τ,divV,Δ,D.ηc,D.ηv,P.g,Fm,FPt)
            F[Num.Vx]   =   Fm.x[:]
            F[Num.Vy]   =   Fm.y[:]
            F[Num.Pt]   =   FPt[:]
            RM          =   norm(F)/length(F)
            if iter == 1
                R0 = RM
            end
            RMrel       =   RM/R0
            @printf("it: %i, ||R|| = %1.4e, ||R||/||R₀|| = %1.4e\n",iter,RM,RMrel)
            (RM < atolM || RM/R0 < rtolM) && break
            end
            # Assemble Coefficients --
            @timeit to "Assembly" begin
            KM      =   Assembly(NC, NV, Δ, D.ηc, D.ηv, VBC, Num)
            end
            # ---
            # Solution of the linear system --
            @timeit to "Solution" begin
            δx      =   - KM \ F
            end
            # ---
            # Update Unknown Variables ---
            D.vx[:,2:end-1]     .+=  δx[Num.Vx]
            D.vy[2:end-1,:]     .+=  δx[Num.Vy]
            D.Pt                .+=  δx[Num.Pt]
            # ---
            # recompute ε and τ for the updated velocity field ---
            ComputeStrainStress2D!(D,ε,τ,divV,Δ)
            # ---
            # Get second invariants ---
            GetSecondInvariants!(ε,τ)
            # ---
            # Update Rheology ---
            UpdateRheo!(ε,Rhe,D,P)
            # ---
        end
        end
        # Get the velocity on the centroids --- 
        # For visualization purposes and Euler advection methods ---
        for i = 1:NC.x
            for j = 1:NC.y
                D.vxc[i,j]  = (D.vx[i,j+1] + D.vx[i+1,j+1])/2
                D.vyc[i,j]  = (D.vy[i+1,j] + D.vy[i+1,j+1])/2
            end
        end
        @. D.vc        = sqrt(D.vxc^2 + D.vyc^2)
        # --------------------------------------------------------------- #
        # Update time stepping ========================================== #
        GetTimeStep!(T,Δ,P,D)
        # --------------------------------------------------------------- #
        # Energy Conservation Equation ================================== #
        # Shear heating ---
        @. D.Hs =   τ.xx .* ε.xx .+ τ.yy .* ε.yy .+ 2.0 .* τ.xyc .* ε.xyc
        # ---
        # Temperature diffusion --- 
        @printf("---Energy Calculation---\n")
        if FD.Method.Diff==:explicit
            ForwardEuler2Dc!(D, P.κ, Δ.x, Δ.y, T.Δ, NC, TBC;
                        ρ₀ = P.ρ₀, cp = P.cp, Qₛ = D.Hs)
        elseif FD.Method.Diff==:implicit
            BackwardEuler2Dc!(D, P.κ, Δ.x, Δ.y, T.Δ, NC, TBC, rhs, K, Num; 
                        ρ₀ = P.ρ₀, cp = P.cp, Qₛ = D.Hs)
        elseif FD.Method.Diff==:CN
            CNA2Dc!(D, P.κ, Δ.x, Δ.y, T.Δ, NC, TBC, rhs, K1, K2, Num; 
                        ρ₀ = P.ρ₀, cp = P.cp, Qₛ = D.Hs)
        elseif FD.Method.Diff==:dc
            D.T0        .=  D.T
            D.T_ex0     .=  D.T_ex
            for iter = 1:niter
                # Evaluate residual
                ComputeResiduals2Dc!(R, D.T, D.T_ex, D.T0, D.T_ex0, ∂2T, 
                        P.κ, TBC, Δ, T.Δ;
                        C = FD.Method.θ, ρ₀ = P.ρ₀, cp = P.cp, Qₛ = D.Hs)
                @printf("||R|| = %1.4e\n", norm(R)/length(R))
                norm(R)/length(R) < ϵ ? break : nothing
                # Assemble linear system
                K  = AssembleMatrix2Dc(P.κ, TBC, Num, NC, Δ, T.Δ;C=FD.Method.θ)
                # Solve for temperature correction: Cholesky factorisation
                Kc = cholesky(K.cscmatrix)
                # Solve for temperature correction: Back substitutions
                δT = -(Kc\R[:])
                # Update temperature
                @. D.T += δT[Num.T]
            end
            D.T_ex[2:end-1,2:end-1]     .=  D.T
        end
        # ---
        # Temperature advection ---
        if FD.Method.Adv==:tracers 
            # Update tracer info from grid ---
            # - Calculate temperature difference new and old time 
            # - Calculate subgrid diffusion on tracers
            # - Interpolate subgrid diffuison to the grid
            # - Calculate remaining temperature difference on the grid 
            # - Interpolate remaining T on tracers
            # - Calculate corrected temperature difference on tracers
        else
            # Advect temperature, Eulerian schemes ---
            # Correct velocity for eulerian advection schemes ---
            @. vx_bg_c  =   -Ini.ε .* x.c2d
            @. vy_bg_c  =   Ini.ε .* y.c2d
            @. vxc_res  =   D.vxc - vx_bg_c 
            @. vyc_res  =   D.vyc - vy_bg_c
            if it == 1
                @. vxco_res    =   D.vxc - vx_bg_c 
                @. vyco_res    =   D.vyc - vy_bg_c 
            end
            if FD.Method.Adv==:upwind
                upwindc2D!(D.T,D.T_ex,vxc_res,vyc_res,NC,T.Δ,Δ.x,Δ.y)
            elseif FD.Method.Adv==:slf
                @show minimum(D.T_exD0), maximum(D.T_exD0)
                @show minimum(D.T_ex0), maximum(D.T_ex0)
                @show minimum(D.T_ex), maximum(D.T_ex)
                # @. D.TD0    =   D.T
                # @. D.T_ex0  =   D.T_ex
                # q = heatmap(x.c./1e3,y.c./1e3,(D.T0'),color=cgrad(:lajolla),
                #     xlabel="x[km]",ylabel="y[km]",
                #     # clims=(-13.5,-12.0),
                #     aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3), 
                #     ylims=(M.ymin/1e3, M.ymax/1e3),colorbar=true,
                #     size = (900,700), titlefontsize = 20,
                #     guidefontsize = 20, tickfontsize = 16,
                #     colorbar_tickfontsize = 14,
                #     colorbar_titlefontsize = 16)
                # display(q)
                slfc2D!(D.T,D.T_ex,D.T_exD0,vxc_res,vyc_res,NC,T.Δ,Δ.x,Δ.y)
                # @. D.T_exD0     =   D.T_ex
                # D.T_ex[2:end-1,2:end-1]  .=   D.T
                # Save diffused field for next SLF step
                copyto!(D.T_exD0, D.T_ex)
                # Replace interior by the advected solution
                D.T_ex[2:end-1,2:end-1] .= D.T
                # q1 = heatmap(x.c./1e3,y.c./1e3,(D.T'),color=cgrad(:lajolla),
                #     xlabel="x[km]",ylabel="y[km]",
                #     # clims=(-13.5,-12.0),
                #     aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3), 
                #     ylims=(M.ymin/1e3, M.ymax/1e3),colorbar=true,
                #     size = (900,700), titlefontsize = 20,
                #     guidefontsize = 20, tickfontsize = 16,
                #     colorbar_tickfontsize = 14,
                #     colorbar_titlefontsize = 16)
                # display(q1)
            elseif FD.Method.Adv==:semilag
                semilagc2D!(D.T,D.T_ex,vxc_res,vyc_res,vxco_res,vyco_res,x,y,T.Δ)
                vxco_res    .=   vxc_res
                vyco_res    .=   vyc_res
            end
        end
        # ---
        # Advection of tracers ---
        @printf("Running on %d thread(s)\n", nthreads())  
        AdvectTracer2D(Ma,nmark,D,x,y,T.Δ,Δ,NC,rkw,rkv)
        # ---
        # Diagnostics =================================================== #
        pp,r    =   Diagnostics!(phbd,halfwidth,nprof,
                        M,T,D,Ini,x,y,ε,
                        θsb,θsb2,Dsb,εf,δTemp,shortening,strain,style,
                        xp,yp,it)
        # Save diagnostics on file ---
        if save_fig == 1
            @printf(dfile,
                "%d %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
                it,
                strain[it],
                shortening[it],
                εf[it],
                δTemp[it],
                Dsb[it],
                θsb[it]*180/π,
                θsb2[it]*180/π)
            flush(dfile)   # immediately write to disk
        end
        # --------------------------------------------------------------- #
        if mod(it,5) == 0 || it == 1 || it == T.itmax
            if save_fig == 1
                pp !== nothing && Plots.frame(animProf, pp)
                r  !== nothing && Plots.frame(anim2D, r)
            else
                pp !== nothing && display(pp)
                r  !== nothing && display(r)
            end
        end
        # D.T0        .=  D.T
        # D.T_ex0     .=  D.T_ex
        # Deform and remesh grid ======================================== #
        # Effectively the new time step --- !!! ---
        # Calculate new xmin, xmax and ymax ---
        M.xmin  =   M.xmin * exp(-Ini.ε * T.Δ)
        M.xmax  =   M.xmax * exp(-Ini.ε * T.Δ)
        M.ymax  =   M.ymax * exp(Ini.ε * T.Δ)
        # Remeshing ---
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
        # --------------------------------------------------------------- #
        # Update boundary velocity ====================================== # 
        IniVelocity!(Ini.V,D,VBC,NV,Δ,M,x,y;Ini.ε)
        # --------------------------------------------------------------- #
        # Update cell centroids and vertices information ================ #
        # if FD.Method.Adv==:tracer
            # Update everything
        # else
            # Update only phase 
        # end
        ϵx = 1e-10 * (M.xmax - M.xmin)
        ϵy = 1e-10 * (M.ymax - M.ymin)
        @. Ma.x = clamp(Ma.x, M.xmin + ϵx, M.xmax - ϵx)
        @. Ma.y = clamp(Ma.y, M.ymin + ϵy, M.ymax - ϵy)
        CountMPC(Ma,nmark,MPC,M,x,y,Δ,NC,NV)
        @timeit to "Tracer Interpolation" begin
        # Interpolate phase from tracers to grid ---
        Markers2Cells(Ma,nmark,MAVG.PC_th,D.p_ex,MAVG.wte_th,D.wte,x,y,Δ,Aparam,phase)
        D.p     .=  D.p_ex[2:end-1,2:end-1]
        # Density ---
        @.  D.ρ     =   P.ρ₀*(1.0 - D.p) + P.ρ₀*D.p
        D.ρ_ex[2:end-1,2:end-1]     .=  D.ρ
        D.ρ_ex[1,:]     .=  D.ρ_ex[2,:]
        D.ρ_ex[end,:]   .=  D.ρ_ex[end-1,:]
        D.ρ_ex[:,1]     .=  D.ρ_ex[:,2]
        D.ρ_ex[:,end]   .=  D.ρ_ex[:,end-1]
        end
        # --------------------------------------------------------------- #
        if mod(it,10) == 0 || it == 1
            @show it
            @show Time[it]
            @show strain[it]
            @show εf[it]
            @show shortening[it]
            @show δTemp[it]
        end
        if shortening[it] >= 35
        # if shortening[it] >= 20
            T.itmax   =   it
            @show it
            @show Time[it]
            @show strain[it]
            @show εf[it]
            @show shortening[it]
            @show δTemp[it]
            @printf("Bulk shortening %g reached maximum value\n",shortening[it])
            break
        else
            T.itmax   =   it
        end
    end # End Time Loop
    end
    if save_fig == 1
        close(dfile)
    end
    # Save Animation ---
    if save_fig == 1
        # Write the frames to a GIF file
        Plots.gif(anim2D, string(path, filename, "_2D.gif"), fps = 10)
        Plots.gif(animProf, string(path, filename, "_profiles.gif"), fps = 10)
        # foreach(rm, filter(startswith(string(path,"00")), readdir(framepath2D,join=true)))
        # foreach(rm, filter(startswith(string(path,"00")), readdir(framepathProf,join=true)))
    end
    t = plot(shortening[1:T.itmax],
                [εf[1:T.itmax]],
                xlims=(0,shortening[T.itmax]),
                label="",ylabel="εf",xlabel="ε [%]",)
    u   =   plot(shortening[1:T.itmax],δTemp[1:T.itmax],
                ylabel="´δT [°C]",xlabel="ε [%]",label="",
                xlims=(0,shortening[T.itmax]),)
    v   =   plot(shortening[1:T.itmax],[θsb[1:T.itmax].*180.0./π θsb2[1:T.itmax].*180.0./π],
                ylabel="θ",xlabel="ε [%]",label="",
                xlims=(0,shortening[T.itmax]),)
    w   =   plot(shortening[1:T.itmax],Dsb[1:T.itmax]./1e3,
                ylabel="D [km]",xlabel="ε [%]",label="",
                xlims=(0,shortening[T.itmax]),)
    if save_fig == 1
        savefig(t,string(path,"/StrainRateMulitply_",FD.Method.Diff,"_",FD.Method.Adv,"_",
                @sprintf("%i",NC.x),"_",@sprintf("%i",NC.y),".png"))
        savefig(u,string(path,"/DeltaTemp_",FD.Method.Diff,"_",FD.Method.Adv,"_",
                @sprintf("%i",NC.x),"_",@sprintf("%i",NC.y),".png"))
        savefig(v,string(path,"/ShearZoneAngle_",FD.Method.Diff,"_",FD.Method.Adv,"_",
                @sprintf("%i",NC.x),"_",@sprintf("%i",NC.y),".png"))
        savefig(w,string(path,"/ShearZoneThickness_",FD.Method.Diff,"_",FD.Method.Adv,"_",
                @sprintf("%i",NC.x),"_",@sprintf("%i",NC.y),".png"))
    else
        display(t)
        display(u)
        display(v)
        display(w)
    end
    display(to)
end
# ======================================================================= #
# =================================== END =============================== # 
# ======================================================================= #

# ======================================================================= #
# Call Main Script ====================================================== #
# ======================================================================= #

ShearHeatingShearBands(:dc,0.5,:slf,:fixed)
# ShearHeatingShearBands()

# ======================================================================= #
# ======================================================================= # 
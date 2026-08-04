using Plots, GeoModBox.HeatEquation.TwoD, ExtendableSparse
using Statistics
using TimerOutputs, LaTeXStrings, Measures
using ExactFieldSolutions

function Gaussian_Diffusion()
to      =   TimerOutput()
Schema  =   ["explicit","implicit","CN","ADI"]
ns          =   size(Schema,1)
nrnxny      =   6
save_fig    =   -1
# Physical Parameters ------------------------------------------------ #
P       = ( 
    L       =   200e3,          #   Length [ m ]
    H       =   200e3,          #   Height [ m ]
    k       =   3,              #   Thermal Conductivity [ W/m/K ]
    cp      =   1000,           #   Specific Heat Capacity [ J/kg/K ]
    ρ       =   3200,           #   Density [ kg/m^3 ]
)
P1      = (
    κ       =   P.k/P.ρ/P.cp,   #   Thermal Diffusivity [ m^2/s ] 
    Tamp    =   500,            #   Temperaturamplitude [K]
    σ       =   20e3,           #   
)
P       =   merge(P,P1)
# -------------------------------------------------------------------- #
# Statistical Parameter ---------------------------------------------- #
St      = (
    ε           =   zeros(size(Schema,1),nrnxny),    
    nxny        =   zeros(size(Schema,1),nrnxny),
    Tmax        =   zeros(size(Schema,1),nrnxny),
    Tmean       =   zeros(size(Schema,1),nrnxny),
)
# -------------------------------------------------------------------- #
# Loop over different discretization schemes ------------------------- #
@timeit to "Discretization Loop" begin
for m = 1:ns
    FDSchema = Schema[m]
    display(FDSchema)
    @timeit to "Resolution Loop" begin
    for l = 1:nrnxny
        @timeit to "Ini" begin
        # Numerical Parameters --------------------------------------- #
        NC  = (
            x       =   l*20,       #   Number of Centroids in x
            y       =   l*20        #   Number of Centroids in y
        )
        Δ   = (
            x       =   P.L/NC.x,   #   Grid spacing in x
            y       =   P.H/NC.y    #   Grid Spacing in y
        )
        display(string("nx = ",NC.x,", ny = ",NC.y))
        # ------------------------------------------------------------ #
        # Animationssettings ----------------------------------------- #
        path        =   string("./examples/DiffusionEquation/2D/Results/")
        anim        =   Plots.Animation(path, String[] )
        filename    =   string("Gaussian_Diffusion_",FDSchema,
                            "_nx_",NC.x,"_ny_",NC.y)
        # ------------------------------------------------------------ #        
        # Grid coordinates ------------------------------------------- #
        x       = (
            c       =   LinRange(-P.L/2+ Δ.x/2.0, P.L/2 - Δ.x/2.0, NC.x),
        )
        y       = (
            c       =   LinRange(-P.H/2 + Δ.y/2.0, P.H/2 - Δ.y/2.0, NC.y),
        )
        # ------------------------------------------------------------ #
        # Time Parameters -------------------------------------------- #
        T       = (
            year        =   365.25*3600*24,     #   Seconds per year
            Δfac        =   1.0,                #   Factor for Explicit Stability Criterion
        )
        T1      = (
            tmax        =   10 * 1e6 * T.year,  #   Maximum Time in [ s ]
            Δ           =   [0.0]            
        )
        T       =   merge(T,T1)
        T.Δ[1]  =   T.Δfac * (1.0 / ( 2.0 * P.κ * ( 1 /Δ.x^2 + 1 / Δ.y^2 )))
        
        nt      =   ceil(Int,T.tmax/T.Δ[1]) + 1     #   Number of Time Steps
        time    =   zeros(nt)
        # ------------------------------------------------------------ #
        # Initial Conditions  ---------------------------------------- #
        D       = (
            T           =   zeros(NC...),
            T0          =   zeros(NC...),
            T_ex        =   zeros(NC.x+2,NC.y+2),
            Tana        =   zeros(NC...),
            RMS         =   zeros(1,nt),
            εT          =   zeros(NC...),
            Tmax        =   zeros(1,nt),
            Tmean       =   zeros(1,nt),
            Tprofile    =   zeros(NC.y,nt),
            Tprofilea   =   zeros(NC.y,nt),           
        )
        # Initial conditions
        AnalyticalSolution2D!(D.T, x.c, y.c, time[1], (T0=P.Tamp,K=P.κ,σ=P.σ))
        @. D.Tana                   =   D.T
        @. D.T0                     =   D.T
        D.T_ex[2:end-1,2:end-1]     .=  D.T
    
        D.Tprofile[:,1]     .=  (D.T[convert(Int,NC.x/2),:] + 
                                    D.T[convert(Int,NC.x/2)+1,:]) / 2
        D.Tprofilea[:,1]    .=  (D.Tana[convert(Int,NC.x/2),:] + 
                                    D.Tana[convert(Int,NC.x/2)+1,:]) / 2
        # Visualize initial condition ---
        # subplot 1 ---
        p = heatmap(x.c ./ 1e3, y.c ./ 1e3, (D.T)', 
                color=:viridis, colorbar=true, aspect_ratio=:equal, 
                xlabel= L"x\ [km]", ylabel= L"z\ [km]", 
                title= L"Temperature\ [K]", 
                size = (900,600), dpi = 300,
                guidefontsize = 12, tickfontsize = 12,
                titlefontsize = 14,
                xlims=(-P.L/2/1e3, P.L/2/1e3), ylims=(-P.H/2/1e3, P.H/2/1e3), 
                clims=(minimum(D.T), maximum(D.T)),layout=(2,2),
                subplot=1)
        contour!(p,x.c./1e3,y.c/1e3,D.T',linewidth=2.0,
                levels=5,linecolor=:black,subplot=1)
        contour!(p,x.c./1e3,y.c/1e3,D.Tana',linewidth=2.0,
                levels=5,linestyle=:dash,linecolor=:yellow,subplot=1)
        annotate!(p,-170, 110,text("a)", 14, :black, :bold),subplot = 1)
        # subplot 2 ---
        heatmap!(p,x.c ./ 1e3, y.c ./ 1e3, D.εT', 
                color=:viridis, colorbar=true, aspect_ratio=:equal, 
                xlabel= L"x\ [km]", ylabel= L"z\ [km]", 
                title= L"Deviation\ [K]", 
                size = (900,600), dpi = 300,
                guidefontsize = 12, tickfontsize = 12,
                titlefontsize = 14,
                xlims=(-P.L/2/1e3, P.L/2/1e3), ylims=(-P.H/2/1e3, P.H/2/1e3),  
                # clims=(-1,1),
                layout=(2,2),
                subplot=2)
        annotate!(p,-170,110,text("b)", 14, :black, :bold),subplot = 2)
        # subplot 3 ---
        plot!(p,D.Tprofile[:,1],y.c./1e3,
                linecolor=:black,linewidth=2.0,
                ylims=(-P.H/2/1e3, P.H/2/1e3),
                xlims=(0,P.Tamp),
                size = (900,600), dpi = 300,
                guidefontsize = 12, tickfontsize = 12,
                bottom_margin = 5mm, left_margin = 5mm,
                xlabel= L"T_{x=0\ km}\ [K]",
                ylabel= L"z\ [km]",
                label="",
                subplot=3)
        plot!(p,D.Tprofilea[:,1],y.c./1e3,
                linestyle=:dash,linecolor=:yellow,linewidth=2.0,
                label="",
                subplot=3)
        annotate!(p,-100,100,text("c)", 14, :black, :bold),subplot = 3)
        # subplot 4 ---
        plot!(p,time[1:end]./T.year./1e6,D.RMS[1:end],
                label="",
                guidefontsize = 12, tickfontsize = 12,
                size = (900,600), dpi = 300,linewidth=2.0,
                xlims = (0,10), ylims = (0,.1),
                bottom_margin = 5mm,
                xlabel= L"t\ [Myrs]",ylabel= L"RMS",
                subplot=4)
        annotate!(p,-1.7,0.1,text("d)", 14, :black, :bold),subplot = 4)
        if save_fig == 0
            display(p)
        end
        # Boundary Conditions ---------------------------------------- #
        BC     = (type    = (W=:Dirichlet, E=:Dirichlet, 
                                N=:Dirichlet, S=:Dirichlet),
                    val     = (W=D.Tana[1,:],E=D.Tana[end,:],
                                N=D.Tana[:,end],S=D.Tana[:,1]))
        # ------------------------------------------------------------ #
        if FDSchema == "implicit"
            # Linear System of Equations ----------------------------- #
            Num     =   (T=reshape(1:NC.x*NC.y, NC.x, NC.y),)
            ndof    =   maximum(Num.T)
            K       =   ExtendableSparseMatrix(ndof,ndof)
            rhs     =   zeros(ndof)
        end
        if FDSchema == "CN"
            # Linear System of Equations ----------------------------- #
            Num     =   (T=reshape(1:NC.x*NC.y, NC.x, NC.y),)
            ndof    =   maximum(Num.T)
            K1      =   ExtendableSparseMatrix(ndof,ndof)
            K2      =   ExtendableSparseMatrix(ndof,ndof)
            rhs     =   zeros(ndof)
        end
        end
        @timeit to "Time Loop" begin
        # Time Loop -------------------------------------------------- #
        for n = 2:nt
            time[n]     =   time[n-1] + T.Δ[1]
            if time[n] > T.tmax 
                T.Δ[1]  =   T.tmax - time[n-1]
                time[n] =   time[n-1] + T.Δ[1]
            end               
            # Exact solution on cell centroids
            AnalyticalSolution2D!(D.Tana, x.c, y.c, time[n], (T0=P.Tamp,K=P.κ,σ=P.σ))
            # Exact solution on cell boundaries
            BoundaryConditions2D!(BC, x.c, y.c, time[n], (T0=P.Tamp,K=P.κ,σ=P.σ)) 
            if FDSchema == "explicit"
                @timeit to "Explicit" begin
                ForwardEuler2Dc!(D, P.κ, Δ.x, Δ.y, T.Δ[1], NC, BC)
                end
            elseif FDSchema == "implicit"
                @timeit to "Implicit" begin
                BackwardEuler2Dc!(D, P.κ, Δ.x, Δ.y, T.Δ[1], NC, BC, rhs, K, Num)
                end
            elseif FDSchema == "CN"
                @timeit to "CN" begin
                CNA2Dc!(D, P.κ, Δ.x, Δ.y, T.Δ[1], NC, BC, rhs, K1, K2, Num)
                end
            elseif FDSchema == "ADI"
                @timeit to "ADI" begin
                ADI2Dc!(D, P.κ, Δ.x, Δ.y, T.Δ[1], NC, BC)
                end
            end
            # Maximum and Mean Temperature with time ---
            D.Tmax[n]   =   maximum(D.T)
            D.Tmean[n]  =   mean(D.T)
            # Vertical Profile along the Center of the Domain ---
            D.Tprofile[:,n]     .=  (D.T[convert(Int,NC.x/2),:] + 
                                        D.T[convert(Int,NC.x/2)+1,:]) / 2
            D.Tprofilea[:,n]    .=  (D.Tana[convert(Int,NC.x/2),:] + 
                                        D.Tana[convert(Int,NC.x/2)+1,:]) / 2
            # Deviation from the Analytical Solution ---
            @. D.εT     =   (D.Tana - D.T)
            # RMS ---
            D.RMS[n]    =   sqrt(sum(D.εT.^2)/(NC.x*NC.y))
            # Plot Solution ---
            if mod(n,2) == 0 || n == nt
                # subplot 1 ---
                p = heatmap(x.c ./ 1e3, y.c ./ 1e3, (D.T)', 
                    color=:viridis, colorbar=true, aspect_ratio=:equal, 
                    xlabel= L"x\ [km]", ylabel= L"z\ [km]", 
                    title= L"Temperature\ [K]", 
                    size = (900,600), dpi = 300,
                    guidefontsize = 12, tickfontsize = 12,
                    titlefontsize = 14,
                    # top_margin = 2mm,
                    # guidefontsize = 22, tickfontsize = 22,
                    xlims=(-P.L/2/1e3, P.L/2/1e3), ylims=(-P.H/2/1e3, P.H/2/1e3), 
                    clims=(minimum(D.T), maximum(D.T)),layout=(2,2),
                    subplot=1)
                contour!(p,x.c./1e3,y.c/1e3,D.T',linewidth=2.0,
                        levels=5,linecolor=:black,subplot=1)
                contour!(p,x.c./1e3,y.c/1e3,D.Tana',linewidth=2.0,
                        levels=5,linestyle=:dash,linecolor=:yellow,subplot=1)
                annotate!(p,-170,110,text("a)", 14, :black, :bold),subplot = 1)
                # subplot 2 ---
                heatmap!(p,x.c ./ 1e3, y.c ./ 1e3, D.εT', 
                    color=:viridis, colorbar=true, aspect_ratio=:equal, 
                    xlabel= L"x\ [km]", ylabel= L"z\ [km]", 
                    title= L"Deviation\ [K]", 
                    size = (900,600), dpi = 300,
                    guidefontsize = 12, tickfontsize = 12,
                    titlefontsize = 14,
                    # top_margin = 2mm,
                    xlims=(-P.L/2/1e3, P.L/2/1e3), ylims=(-P.H/2/1e3, P.H/2/1e3),  
                    # clims=(-1,1),
                    layout=(2,2),
                    subplot=2)
                annotate!(p,-170,110,text("b)", 14, :black, :bold),subplot = 2)
                # subplot 3 ---
                plot!(p,D.Tprofile[:,n],y.c./1e3,
                        linecolor=:black,linewidth=2.0,
                        ylims=(-P.H/2/1e3, P.H/2/1e3),
                        xlims=(0,P.Tamp),
                        size = (900,600), dpi = 300,
                        guidefontsize = 12, tickfontsize = 12,
                        bottom_margin = 5mm, left_margin = 5mm,
                        xlabel= L"T_{x=0\ km}\ [K]",
                        ylabel= L"z\ [km]",
                        label="",
                        subplot=3)
                plot!(p,D.Tprofilea[:,n],y.c./1e3,
                        linestyle=:dash,linecolor=:yellow,linewidth=2.0,
                        label="",
                        subplot=3)
                annotate!(p,-100,100,text("c)", 14, :black, :bold),subplot = 3)
                # subplot 4 ---
                plot!(p,time[1:n]./T.year./1e6,D.RMS[1:n],
                    label="",linewidth=2.0,
                    size = (900,600), dpi = 300,
                    guidefontsize = 12, tickfontsize = 12,
                    bottom_margin = 5mm,
                    xlims = (0,10), ylims = (0,.1),
                    xlabel= L"t\ [Myrs]",ylabel= L"RMS",
                    subplot=4)
                annotate!(p,-1.7,0.1,text("d)", 14, :black, :bold),subplot = 4)
                if save_fig == 1
                    Plots.frame(anim)
                elseif save_fig == 0
                    display(p)                        
                end
            end
            # End Time Loop ---
        end        
        display("Time loop finished ...")
        display("-> Use new grid size...")
        end
        # Save Animation ---
        if save_fig == 1
            # Write the frames to a GIF file
            Plots.gif(anim, string( path, filename, ".gif" ), fps = 15)
        elseif save_fig == 0
            display(plot(p))
        end
        foreach(rm, filter(startswith(string(path,"00")), readdir(path,join=true)))
        # ------------------------------------------------------------ #
        # Statistical Values for Each Scheme and Resolution ---
        St.ε[m,l]       =   maximum(D.RMS[:])
        St.nxny[m,l]    =   1/NC.x/NC.y
        St.Tmax[m,l]    =   D.Tmax[nt]
        St.Tmean[m,l]   =   D.Tmean[nt]
        # ------------------------------------------------------------ #
    end
    end
end
end
# Visualize Statistical Values --------------------------------------- #
q   =   plot(0,0,layout=(1,3),
            dpi=300) 
for m = 1:ns
    plot!(q,St.nxny[m,:],St.ε[m,:],
                marker=:circle,markersize=4,
                legend = :topleft,
                label=Schema[m],
                xaxis=:log,yaxis=:log,
                markerstrokewidth=0.0,
                xlims=(3e-5,5e-3),
                ylims=(1e-2,1e1),
                xlabel= L"\frac{1}{nx \cdot ny}",ylabel= L"ε_{T}",
                subplot=1)
    plot!(q,St.nxny[m,:],St.Tmax[m,:],
                marker=:circle,markersize=4,label="",
                xaxis=:log,
                xlims=(3e-5,5e-3),
                ylims=(86,100),
                markerstrokewidth=0.0,
                xlabel=L"\frac{1}{nx \cdot ny}",ylabel= L"T_{max}",
                subplot=2)
    plot!(q,St.nxny[m,:],St.Tmean[m,:],
                marker=:circle,markersize=4,label="",
                xaxis=:log,
                xlims=(3e-5,5e-3),
                ylims=(9.97,10.01),
                markerstrokewidth=0.0,
                xlabel=L"\frac{1}{nx \cdot ny}",ylabel= L"⟨\ T\ ⟩",
                subplot=3)
end
annotate!(
    q, 3.5e-6, 10.0,
    text("a)", 10, :black, :bold, :left),
    subplot = 1,
)
annotate!(
    q, 3.5e-6, 100.0,
    text("b)", 10, :black, :bold, :left),
    subplot = 2,
)
annotate!(
    q, 3.5e-6, 10.01,
    text("c)", 10, :black, :bold, :left),
    subplot = 3,
)
display(q)
# --------------------------------------------------------------------- #
# Save Final Figure --------------------------------------------------- #
if save_fig == -1 || save_fig == 1
    savefig(q,"./examples/DiffusionEquation/2D/Results/Gaussian_ResTest.png")
end
# --------------------------------------------------------------------- #
display(to)
end

Gaussian_Diffusion()


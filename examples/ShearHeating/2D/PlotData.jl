using DelimitedFiles
using Plots
using LaTeXStrings
using FileIO, ImageIO

# ============================================================
# Helper function 
# ============================================================
function field_plot(fieldnr, basepath, Diff, θ, Adv, NCx, NCy, avg_p, avg_v, style; title="")

    folder = joinpath(
        basepath,
        string(Diff, "_", θ, "_", Adv, "_", NCx, "_", NCy, "_",
        avg_p, "_", avg_v, "_", style),
        "frames_2D/",
    )

    file = joinpath(folder, string(fieldnr, ".png"))
    # file = string(folder,fieldnr,".png")
    if !isfile(file)
        error("Field image not found: $(file)")
    end

    @show file
    img = load(file)
    # @show img

    return plot(
        img,
        axis = nothing,
        ticks = nothing,
        border = :none,
        legend = false,
        title = title,
        aspect_ratio = :equal,
    )
end


function read_diagnostics(file)

    data = readdlm(file, comments=true, comment_char='#')

    return (
        it         = data[:,1],
        strain     = data[:,2],
        shortening = data[:,3],
        εf         = data[:,4],
        δT         = data[:,5],
        Dsb        = data[:,6],
        θsb        = data[:,7],
        θc         = data[:,8],
    )

end

function diagnostics_file(basepath, Diff, θ, Adv, NCx, NCy, avg_p, avg_v, style)

    folder = joinpath(
        basepath,
        string(Diff, "_", θ, "_", Adv, "_", NCx, "_", NCy, "_",
        avg_p, "_", avg_v, "_", style)
    )

    return joinpath(
        folder,
        string(
            "Diagnostics_",
            Diff, "_",
            θ, "_",
            Adv, "_",
            NCx, "_",
            NCy, "_.txt"
        )
    )

end


# ============================================================
# Settings
# ============================================================

basepath = "./examples/ShearHeating/2D/Results"

Diff = :dc

# θlist = [0.0, 0.5, 1.0]
θlist = [0.5]

Advlist = [
    # :upwind,
    # :slf,
    :semilag,
]

style = :fixed

avg_p = :geometric
avg_v = :geometric

NCx, NCy = 200, 100

outpath = joinpath(basepath, "DiagnosticsPlots")
isdir(outpath) || mkpath(outpath)

# ============================================================
# Read Duretz et al. Figure 2 (digitized)
# ============================================================

duretz = readdlm("./examples/ShearHeating/2D/Duretz_Figure2_digitized.txt",
                 comments=true,
                 comment_char='#')

short_D =   duretz[:,1]
D_D     =   duretz[:,2]
θ_D     =   duretz[:,3]
δT_D    =   duretz[:,4]
εf_D    =   duretz[:,5]

# ============================================================
# Plot styles
# ============================================================

colors = Dict(
    :upwind  => :royalblue,
    :slf     => :black,
    :semilag => :firebrick,
)

linestyles = Dict(
    0.0 => :solid,
    0.5 => :dash,
    1.0 => :dot,
)

linewidth = 2.5

# ============================================================
# Create empty plots
# ============================================================

p1 = plot(
    xlabel = "Bulk shortening [%]",
    ylabel = L"\varepsilon_f",
    title = "Strain-rate amplification",
    legend = :bottomright,
    framestyle = :box,
)

p2 = plot(
    xlabel = "Bulk shortening [%]",
    ylabel = L"\delta T\ [^\circ C]",
    title = "Temperature increase",
    legend = false,
    framestyle = :box,
)

p3 = plot(
    xlabel = "Bulk shortening [%]",
    ylabel = L"\theta_{SB}\ [^\circ]",
    title = "Shear-band angle",
    legend = false,
    framestyle = :box,
)

p4 = plot(
    xlabel = "Bulk shortening [%]",
    ylabel = L"D_{SB}\ [km]",
    title = "Shear-band thickness",
    legend = false,
    framestyle = :box,
)

# ============================================================
# ----- Dummy legend -----
# ============================================================

# blank line

plot!(p1, [NaN], [NaN],
    linecolor=:white,
    label="Advection")

plot!(p1, [NaN], [NaN],
    color=colors[:upwind],
    linewidth=linewidth,
    label="Upwind")

plot!(p1, [NaN], [NaN],
    color=colors[:slf],
    linewidth=linewidth,
    label="Staggered Leapfrog")

plot!(p1, [NaN], [NaN],
    color=colors[:semilag],
    linewidth=linewidth,
    label="Semi-Lagrangian")

# blank line

plot!(p1, [NaN], [NaN],
    linecolor=:white,
    label="")

plot!(p1, [NaN], [NaN],
    linecolor=:white,
    label="Diffusion")

plot!(p1, [NaN], [NaN],
    color=:black,
    linestyle=:solid,
    linewidth=linewidth,
    label="θ = 0.0")

plot!(p1, [NaN], [NaN],
    color=:black,
    linestyle=:dash,
    linewidth=linewidth,
    label="θ = 0.5")

plot!(p1, [NaN], [NaN],
    color=:black,
    linestyle=:dot,
    linewidth=linewidth,
    label="θ = 1.0")

# blank line
plot!(p1, [NaN], [NaN],
    linecolor=:white,
    label="")

scatter!(p1, [NaN], [NaN],
    color=:black,
    marker=:circle,
    markersize=3,
    markerstrokewidth=0,
    label="Duretz et al. (2014)")

# ============================================================
# Plot all simulations
# ============================================================

θfield = θlist[1]
Advfield = Advlist[1]

f1 = field_plot("000007", basepath, Diff, θfield, Advfield,
                NCx, NCy, avg_p, avg_v, style; title="ε = 5%")

f2 = field_plot("000020", basepath, Diff, θfield, Advfield,
                NCx, NCy, avg_p, avg_v, style; title="ε = 15%")

f3 = field_plot("000037", basepath, Diff, θfield, Advfield,
                NCx, NCy, avg_p, avg_v, style; title="ε = 25%")

for θ in θlist

    for Adv in Advlist

        # @show f1
        file = diagnostics_file(basepath, Diff, θ, Adv, NCx, NCy, avg_p, avg_v, style)

        if !isfile(file)
            @warn "File not found" file
            continue
        end

        d = read_diagnostics(file)

        plot!(
            p1,
            d.shortening,
            d.εf,
            color = colors[Adv],
            linestyle = linestyles[θ],
            linewidth = linewidth,
            label = false,
        )

        plot!(
            p2,
            d.shortening,
            d.δT,
            color = colors[Adv],
            linestyle = linestyles[θ],
            linewidth = linewidth,
            label = false,
        )

        plot!(
            p3,
            d.shortening,
            d.θsb,
            color = colors[Adv],
            linestyle = linestyles[θ],
            linewidth = linewidth,
            label = false,
        )

        plot!(
            p4,
            d.shortening,
            d.Dsb ./ 1000,
            color = colors[Adv],
            linestyle = linestyles[θ],
            linewidth = linewidth,
            label = false,
        )

    end

end

# ============================================================
# Duretz et al. (2016)
# ============================================================

for (plt,x,y) in (
    (p1, short_D, εf_D),
    (p2, short_D, δT_D),
    (p3, short_D, θ_D),
    (p4, short_D, D_D),
)
    # plot!(plt, x, y,
    #     color=:black,
    #     linewidth=1,
    #     label=false)

    scatter!(plt, x, y,
        color=:black,
        marker=:circle,
        markersize=4,
        markerstrokewidth=0,
        label=false)
end

# ============================================================
# Axis limits
# ============================================================

ylims!(p1,(0,15))
xlims!(p1,(0,30))
ylims!(p2,(0,170))
xlims!(p2,(0,30))
ylims!(p3,(45,75))
xlims!(p3,(0,30))
ylims!(p4,(0,7))
xlims!(p4,(0,30))

# ============================================================
# Combined figure
# ============================================================

fig = plot(
    p1,p2,
    p3,p4,
    layout=(2,2),
    size=(1200,900)
)
fig2 = plot(
    f1, p1,
    f2, p2,
    f3, p4,
    layout = (3,2),
    size = (1200,1400),
)

display(fig)

# savefig(fig, joinpath(outpath,"DiagnosticsComparison.png"))

fig2 = plot(
    f1, p1,
    f2, p2,
    f3, p4,
    layout = (3,2),
    size = (1200,1400),
)

display(fig2)

# savefig(fig2, joinpath(outpath, "FinalBenchmarkFigure.png"))
using DelimitedFiles
using Plots
using LaTeXStrings
using FileIO, ImageIO, Measures

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
    if !isfile(file)
        error("Field image not found: $(file)")
    end

    @show file
    img = load(file)

    return plot(
        img,
        axis = nothing,
        ticks = nothing,
        border = :none,
        legend = false,
        title = title,
        titlefontsize = 22,
        aspect_ratio = :equal,
        dpi = 300,
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

function add_panel!(p, label, labelsize; dx=0.16, dy=0.0)

    xmin, xmax = xlims(p)
    ymin, ymax = ylims(p)

    x = xmin - dx*(xmax-xmin)
    y = ymax - dy*(ymax-ymin)

    annotate!(p, x, y,
        text(label, labelsize, :black, :left))

end

# ============================================================
# Settings
# ============================================================

basepath = "./examples/ShearHeating/2D/Results"

Diff = :dc

θlist = [0.0, 0.5, 1.0]

Advlist = [
    :upwind,
    :semilag,
    :tracers,
]

style = :fixed

avg_p = :arithmetic
avg_v = :arithmetic

NCx, NCy = 200, 100

outpath = joinpath(basepath, "DiagnosticsPlots")
isdir(outpath) || mkpath(outpath)

θfield      =   θlist[2]
Advfield    =   Advlist[3]

# ============================================================
# Plot styles
# ============================================================

colors = Dict(
    :upwind  => :royalblue,
    :semilag => :black,
    :tracers => :firebrick,
)

linestyles = Dict(
    0.0 => :solid,
    0.5 => :dash,
    1.0 => :dot,
)

linewidth = 3.0

# ============================================================
# Font sizes
# ============================================================

guidefontsize  = 20
tickfontsize   = 16
legendfontsize = 14

# ============================================================
# Create empty plots
# ============================================================

p1 = plot(
    ylabel          = L"\varepsilon_f",
    title           = "",
    legend          = :bottomright,
    framestyle      = :box,
    guidefontsize   = guidefontsize,
    tickfontsize    = tickfontsize,
    legendfontsize  = legendfontsize,
    xlabel          = "",
    xformatter      = _ -> ""
)

p2 = plot(
    ylabel          = L"\delta T\ [^\circ C]",
    title           = "",
    legend          = false,
    framestyle      = :box,
    guidefontsize   = guidefontsize,
    tickfontsize    = tickfontsize,
    legendfontsize  = legendfontsize,
    xlabel          = "",
    xformatter      = _ -> ""
)

p3 = plot(
    xlabel          = "Bulk shortening [%]",
    ylabel          = L"\theta_{SB}\ [^\circ]",
    title           = "",
    legend          = false,
    framestyle      = :box,
    guidefontsize   = guidefontsize,
    tickfontsize    = tickfontsize,
    legendfontsize  = legendfontsize,
)

p4 = plot(
    xlabel          = "Bulk shortening [%]",
    ylabel          = L"D_{SB}\ [km]",
    title           = "",
    legend          = false,
    framestyle      = :box,
    guidefontsize   = guidefontsize,
    tickfontsize    = tickfontsize,
    legendfontsize  = legendfontsize,
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
    color=colors[:tracers],
    linewidth=linewidth,
    label="Tracers")

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

# ============================================================
# Plot all simulations
# ============================================================

f1 = field_plot("000007", basepath, Diff, θfield, Advfield,
                NCx, NCy, avg_p, avg_v, style; title="ε ≈ 5%")

f2 = field_plot("000020", basepath, Diff, θfield, Advfield,
                NCx, NCy, avg_p, avg_v, style; title="ε ≈ 15%")

f3 = field_plot("000037", basepath, Diff, θfield, Advfield,
                NCx, NCy, avg_p, avg_v, style; title="ε ≈ 25%")

for θ in θlist

    for Adv in Advlist

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
            right_margin = 5mm,
        )

        plot!(
            p2,
            d.shortening,
            d.δT,
            color = colors[Adv],
            linestyle = linestyles[θ],
            linewidth = linewidth,
            label = false,
            right_margin = 5mm,
        )

        plot!(
            p3,
            d.shortening,
            d.θsb,
            color = colors[Adv],
            linestyle = linestyles[θ],
            linewidth = linewidth,
            label = false,
            right_margin = 5mm,
        )

        plot!(
            p4,
            d.shortening,
            d.Dsb ./ 1000,
            color = colors[Adv],
            linestyle = linestyles[θ],
            linewidth = linewidth,
            label = false,
            right_margin = 5mm,
        )

    end

end

# ============================================================
# Axis limits
# ============================================================

ylims!(p1,(0,20))
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
    size=(1600,1600)
)

display(fig)

savefig(fig, joinpath(outpath,string("DiagnosticsComparison_",
            style,"_",avg_p,"_",avg_v,".png")))

# ============================================================
# Panel labels
# ============================================================

labelsize = 18

add_panel!(f1, "a)",labelsize;dx=-0.1,dy=1.0)
add_panel!(f2, "b)",labelsize;dx=-0.1,dy=1.0)
add_panel!(f3, "c)",labelsize;dx=-0.1,dy=1.0)
add_panel!(p1, "d)",labelsize)
add_panel!(p2, "e)",labelsize)
add_panel!(p4, "f)",labelsize)

fig2 = plot(
    f1, p1,
    f2, p2,
    f3, p4,
    layout = (3,2),
    size=(1600,1650),
    dpi = 300, 
)

display(fig2)

savefig(fig2, joinpath(outpath,string("FinalBenchmarkFigure_",
                style,"_",avg_p,"_",avg_v,".png")))
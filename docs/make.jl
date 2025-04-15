using Documenter
using GeoModBox

GMB_root_dir = dirname(@__DIR__)

license = read(joinpath(GMB_root_dir, "LICENSE"), String)
write(joinpath(@__DIR__, "src", "man", "license.md"), license)

makedocs(
    sitename = "GeoModBox.jl",
    format = Documenter.HTML(),
    warnonly = [:missing_docs],
    modules = [GeoModBox],
    pages = [
        "Home" => "index.md",
        "Governing Equation" => Any[
            # "Solution" => "man/GESolution.md",
            "Heat Diffusion Equation" => Any[
                "General" => "man/DiffMain.md",
                "1D" => "man/DiffOneD.md",
                "2D" => "man/DiffTwoD.md",
            ],
            # "Advection Equation" => Any[
            #     "General" => "man/AdvMain.md",
            #     "1D" => "man/AdvOneD.md",
            #     "2D" => "man/AdvTwoD.md",
            # ],
            "Momentum Equation" => Any[
                "General" => "man/MomentumMain.md",
                "1D" => "man/MomentumOneD.md",
                "2D" => "man/MomentumTwoD.md",
            ],
        ],
        # "Examples" => "man/Examples.md",
        "List of functions" => "man/listoffunctions.md",
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
withenv("GITHUB_REPOSITORY" => "GeoSci-FFM/GeoModBox.jl") do
    deploydocs(
        # repo = "https://lukasfuchs.github.io/GeoModBox.jl",
        repo = "github.com/GeoSci-FFM/GeoModBox.jl",
        devbranch = "main",
        push_preview = true
    )
end

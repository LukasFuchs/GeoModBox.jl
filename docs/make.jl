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
            "Solution" => Any[
               "General" => "man/GESolution.md", # checked!
            #   "Initial Condition" => "man/Ini.md",
            ],
            "Heat Diffusion Equation" => Any[
                "General" => "man/DiffMain.md", # checked!
                "1D" => "man/DiffOneD.md", # checked!
                "2D" => "man/DiffTwoD.md", 
            ],
            # "Advection Equation" => Any[
            #     "General" => "man/AdvMain.md",
            #     "1D" => "man/AdvOneD.md",
            #     "2D" => "man/AdvTwoD.md",
            # ],
            "Momentum Equation" => Any[
                "General" => "man/MomentumMain.md", # Checked! 
                "1D" => "man/MomentumOneD.md", # Checked! 
                "2D" => "man/MomentumTwoD.md", # Checked!
            ],
        ],
        "Examples" => Any[
              "Overview" => "man/Examples.md",
              "Oceanic Geotherm" => "man/OceanicGeotherm.md",
        ],
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

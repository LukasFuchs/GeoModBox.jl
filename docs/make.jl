using Documenter
using GeoModBox


GMB_root_dir = dirname(@__DIR__)

license = read(joinpath(GBM_root_dir, "LICENSE.md"), String)
write(joinpath(@__DIR__, "src", "man", "license.md"), license)


makedocs(
    sitename = "GeoModBox",
    format = Documenter.HTML(),
    modules = [GeoModBox],
    pages = [
        "Home" => "index.md",
        "Equation" => Any[
            "Energy Conservation" => Any[
                "Diffusion" => "man/DiffMain.md",
            ]
        ],
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#

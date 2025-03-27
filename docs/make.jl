using Documenter
using GeoModBox


GMB_root_dir = dirname(@__DIR__)

license = read(joinpath(GMB_root_dir, "LICENSE"), String)
write(joinpath(@__DIR__, "src", "man", "license.md"), license)


makedocs(
    sitename = "GeoModBox",
    format = Documenter.HTML(),
    warnonly = [:missing_docs],
    modules = [GeoModBox],
    pages = [
        "Home" => "index.md",
        "Equation" => Any[
            "Energy Conservation" => Any[
                "Diffusion" => "man/DiffMain.md",
            ]
        ],
        "List of functions" => "man/listoffunctions.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "https://lukasfuchs.github.io/GeoModBox.jl/",
    devbranch = "main",
    push_preview = true
)

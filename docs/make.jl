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
            "Solution" => "man/GESolution.md",
            "Heat Diffusion Equation" => Any[
                "General" => "man/DiffMain.md",
                "1D" => "man/DiffOneD.md",
                # "2D" => "man/DiffTwoD.md",
            ],
            # "Advection Equation" => Any[
            #     "General" => "man/AdvMain.md",
            #     "1D" => "man/AdvOneD.md",
            #     "2D" => "man/AdvTwoD.md",
            # ],
            # "Momentum Equation" => Any[
            #     "General" => "man/MomMain.md",
            #     "1D" => "man/MomOneD.md",
            #     "2D" => "man/MomTwoD.md",
            # ],
        ],
        # "Examples and Benchmarks" => any[
        #       "Advection" => Any[
        #           
        #       ],
        #       "Diffusion" => Any[
        #           
        #       ],
        #       "Mixed Heated Convection" => Any[
        #           
        #       ],
        #       "Stokes Equation" => Any[
        #           
        #       ],
        # ],
        # "Exercises" => any[
        #       "01 - Euler Advection" => "man/*md.",
        #       "02 - Heat Diffusion, Explicit (1D)" => "man/*md.",
        #       "03 - Heat Diffusion, Implicit (1D)" => "man/*md.",
        #       "04 - Stationary Heat Diffusion (2D)" => "man/*md.",
        #       "05 - Time-Dependent Heat Diffusion (2D)" => "man/*md.",
        #       "06 - Advection (1D)" => "man/*md.",
        #       "07 - Energy Equation (2D)" => "man/*md.",
        #       "08 - Stokes (1D)" => "man/*md.",
        #       "09 - Falling Block (2D), Instantaneous" => "man/*md.",
        #       "10 - Falling Block (2D), Time-Dependent" => "man/*md.",
        # ],
        "List of functions" => "man/listoffunctions.md",
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
withenv("GITHUB_REPOSITORY" => "LukasFuchs/GeoModBox.jl") do
    deploydocs(
        # repo = "https://lukasfuchs.github.io/GeoModBox.jl",
        repo = "github.com/LukasFuchs/GeoModBox.jl",
        devbranch = "main",
        push_preview = true
    )
end

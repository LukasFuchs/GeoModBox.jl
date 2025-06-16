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
            "Advection Equation" => Any[
                 "General" => "man/AdvectMain.md", # Checked! 
                 "1D" => "man/AdvOneD.md", # Checked! 
            #     "2D" => "man/AdvTwoD.md",
            ],
            "Momentum Equation" => Any[
                "General" => "man/MomentumMain.md", # Checked! 
                "1D" => "man/MomentumOneD.md", # Checked! 
                "2D" => "man/MomentumTwoD.md", # Checked!
            ],
        ],
        "Examples" => Any[
            # "Overview" => "man/Examples.md",
            "Diffusion Equation" => Any[
                "Oceanic Geotherm (1D)" => "man/examples/OceanicGeotherm.md",
                "Continental Geotherm (1D)" => "man/examples/ContinentalGeotherm.md",
                "Gaussian Diffusion (1D)" => "man/examples/GaussianDiffusion1D.md",
                "Backward Euler (2D)" => "man/examples/BackwardEuler_DC.md",
                "Forward Euler (2D)" => "man/examples/ForwardEuler_DC.md",
                "Gaussian Diffusion (2D)" => "man/examples/GaussianDiffusion2D.md",
                "Poisson Problem (2D)" => "man/examples/PoissonRestest.md",
                "Poisson Problem; variable k (2D)" => "man/examples/PoissonVariablek.md",
            ],
            # "Advection Equtaion" => Any[
            #     "Advection (2D)" => "man/examples/Advection2D.md",
            #     "Advection Resolution Test (2D)" => "man/examples/AdvectionRestest2D.md",
            # ],
            # "Stokes Equation" => Any[
            #     "Channel Flow (1D)" => "man/examples/ChannelFlow1D.md",
            #     "Falling Block Benchmark" => "man/examples/FallingBlockBenchmark.md",
            #     "Falling Block" => "man/examples/FallingBlockDC.md",
            #     "Rayleigh Taylor Instability" => "man/examples/RTI.md",
            #     "Viscous Inclusion" => "man/examples/ViscousInclusion.md",
            # ],
            # "Mixed Thermal Convection" => Any[
            #     "Bottom Heated" => "man/examples/BottomHeatedConvection.md",
            #     "Internally Heated" => "man/examples/InternallyHeatedConvection.md",
            #     "Mixed Heated" => "man/examples/MixedHeatedConvection.md",
            # ],
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

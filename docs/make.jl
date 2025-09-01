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
        "Home" => "index.md", # Checked! 
        "Governing Equation" => Any[ # Checked! 
            "Solution" => Any[
               "General" => "man/GESolution.md", # checked!
               "Initial Condition" => "man/Ini.md", # Checked!
            ],
            "Heat Diffusion Equation" => Any[
                "General" => "man/DiffMain.md", # checked!
                "1D" => "man/DiffOneD.md", # checked!
                "2D" => "man/DiffTwoD.md", # checked! 
            ],
            "Advection Equation" => Any[
                 "General" => "man/AdvectMain.md", # Checked! 
                 "1D" => "man/AdvOneD.md", # Checked! 
                 "2D" => "man/AdvTwoD.md", # Checked!
            ],
            "Momentum Equation" => Any[
                "General" => "man/MomentumMain.md", # Checked! 
                "1D" => "man/MomentumOneD.md", # Checked! 
                "2D" => "man/MomentumTwoD.md", # Checked!
            ],
        ],
        "Exercises" => Any[
            "01 - Euler Advection" => "man/exercises/01_Euler_Advection.md",
            "02 - Heat Diffusion (explicit)" => "man/exercises/02_1D_Heat_explicit.md",
        ],
        "Examples" => Any[
            "General" => "man/Examples.md",
            "Diffusion Equation" => Any[
                "Oceanic Geotherm (1D)" => "man/examples/OceanicGeotherm.md", # Checked! 
                "Continental Geotherm (1D)" => "man/examples/ContinentalGeotherm.md", # Checked! 
                "Gaussian Diffusion (1D)" => "man/examples/GaussianDiffusion1D.md", # Checked!  
                "Backward Euler (2D)" => "man/examples/BackwardEuler_DC.md", # Checked!
                "Forward Euler (2D)" => "man/examples/ForwardEuler_DC.md", # Checked! 
                "Gaussian Diffusion (2D)" => "man/examples/GaussianDiffusion2D.md", # Checked!
                "Poisson Problem (2D)" => "man/examples/PoissonRestest.md", # Checked!
                "Poisson Problem; variable k (2D)" => "man/examples/PoissonVariablek.md", # Checked!
            ],
            "Advection Equtaion" => Any[
                "Advection (2D)" => "man/examples/Advection2D.md", # Checked! 
                "Advection Resolution Test (2D)" => "man/examples/AdvectionRestest2D.md", # Checked! 
            ],
            "Stokes Equation" => Any[
                "Channel Flow (1D)" => "man/examples/ChannelFlow1D.md", # Checked! 
                "Falling Block Benchmark" => "man/examples/FallingBlockBenchmark.md", # Checked! 
                "Falling Block" => "man/examples/FallingBlockDC.md", # Checked! 
                "Rayleigh Taylor Instability (RTI)" => "man/examples/RTI.md", # Checked!
                "RTI - Growth Rate" => "man/examples/RTI_growth_rate.md", # Checked!
                "Viscous Inclusion" => "man/examples/ViscousInclusion.md",
            ],
            "Mixed Thermal Convection" => Any[
                "Bottom Heated" => "man/examples/BottomHeatedConvection.md",
                "Internally Heated" => "man/examples/InternallyHeatedConvection.md",
                "Mixed Heated" => "man/examples/MixedHeatedConvection.md",
            ],
        ],
        "List of functions" => "man/listoffunctions.md",
        "License" => "man/license.md",
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

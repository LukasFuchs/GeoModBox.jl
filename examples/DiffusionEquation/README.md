# Energy Conservation Equation

This directory contains examples to solve the *diffusive part* of the *temperature equation* in 1-D and 2-D. The examples include

- the determination of an [oceanic](./1D/OceanicGeotherm_1D.jl) and [continental](./1D/ContinentalGeotherm_1D.jl) 1-D geotherm profile, 
- [a comparison of the different **FD**-schemes applied on a 1-D gaussian temperature anomaly](./1D/Heat_1D_discretization.jl), 
- [a backward](./2D/BackwardEuler.jl) and [forward](./2D/ForwardEuler.jl) Euler solution using the defection correction method
- [a 2-D resolution test for each **FD**-scheme using a gaussian temperature anomaly](./2D/Gaussian_Diffusion.jl), 
- [a 2-D poisson problem assuming a variable conductivity](./2D/Poisson_variable_k.jl), and
- [a resolution test for a 2-D poisson problem](./2D/Poisson_ResTest.jl).

For more details on the different discretization methods and each example, please, see the [documentation](https://lukasfuchs.github.io/GeoModBox.jl/).
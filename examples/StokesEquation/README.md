# Momentum Conservation Equation

This directory contains examples to solve the *stokes equation* in 1-D and 2-D. The examples include: 

The examples for a pure stokes equation problem include: 
- [a 1-D channel flow problem for constant and depth-dependent viscosity](./1D/ChannelFlow_1D.jl)

- [a 2-D falling block benchmark](./2D/FallingBlockBenchmark.jl)

- [a 2-D falling block example with constant viscosity using the defect correction method](./2D/FallingBlockConstEta_DC.jl)

- [a 2-D falling block example with variable viscosity using the defect correction method](./2D/FallingBlockVarEta_DC.jl)

- [a 2-D viscous inclusion problem](./2D/ViscousInclusion.jl), and

- [a 2-D Rayleigh-Taylor instability benchmark](./2D/RTI.jl)

For more details on how to solve the *stokes equation* using different methods, please see the [documentation](https://lukasfuchs.github.io/GeoModBox.jl/).
# Mixed Heated Convection

This directory contains examples of three types of mixed heated thermal convection, each driven by distinct heating mechanisms:

- [Purely internally heated](InternallyHeated.jl)
- [Purely basally heated](BottomHeated.jl)
- [Internally and basally heated](MixedHeated.jl)

The default basal *Rayleigh number* $Ra_B$ is set to $10^4$ for all examples; however, higher values of $Ra_B$ (tested up to $Ra_B = 10^6$) can be applied.

All simulations use non-dimensional quantities. The governing equations are solved as follows:
- The momentum conservation equation is solved using the defect correction method.
- The energy conservation equation is solved using the Crank–Nicolson approach.
- Temperature advection is performed using a semi-Lagrangian scheme.

Each model runs for up to 6000 iterations or until a steady state is reached.

For further details on the numerical methods used for each governing equation and for additional information on each example, please refer to the [GeoModBox.jl documentation](https://geosci-ffm.github.io/GeoModBox.jl/).

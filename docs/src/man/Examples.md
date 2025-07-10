# Examples and Benchmarks

This section provides various one- and two-dimensional examples and benchmark problems for each of the governing equations. The examples demonstrate how to implement different numerical solvers, apply scaling, and evaluate the advantages and limitations of various finite difference schemes.

By clicking on the title of each page, you will be directed to the corresponding Julia script in the [examples directory](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples).

> **Note:** In `GeoModBox.jl`, thermal and kinematic boundary conditions are explicitly implemented within the solvers. The absolute values at *ghost nodes* are computed based on the values provided in the boundary condition tuple `BC`. Each tuple specifies the `type` (either Dirichlet or Neumann) and the corresponding `val`ue at each boundary.  
> For constant velocity boundary conditions, additional values must be defined in `val` for the boundary nodes (e.g., `BC.val.vxW`, `BC.val.vxE`). These additional values are required to directly solve the momentum equation and update the right-hand side of the linear system. Furthermore, these values must be assigned to the initial boundary nodes of the respective velocity fields.  
> For more details on the implementation of constant velocity boundaries, refer to the documentation of the [viscous inclusion](./examples/ViscousInclusion.md) example or the [initial velocity setup](Ini.md).

> **Note:** By default, the results of time-dependent examples in `GeoModBox.jl` are stored as GIF animations. To visualize solutions at specific time steps without generating a GIF, set the parameter `save_fig = 0`. In this case, individual plots are not saved, so caution is advised when running problems that require multiple time step iterations.

> **Note:** Some examples use *named tuples* to define constants and parameters. Alternatively, *mutable structures* can be used—particularly useful when parameters need to be modified after initialization (e.g., for scaling purposes). A full transition from *named tuples* to *mutable structures* is planned for future versions of `GeoModBox.jl`.

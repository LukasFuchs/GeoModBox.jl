# Examples and Benchmarks

In the following, one can find different examples and benchmarks in 1-D and 2-D for each of the governing equations. The examples highlight how to implement the different solvers, how to utilize the scaling, and advantages and disadvantages of each finite difference scheme. 

By clicking on the title of each page one is directly directed to the julia example file within the [example directory](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/examples)

> **Note:** Within ```GeoModBox.jl``` the thermal and kinematic boundary conditions are explicitly implemented in the solvers, where the absolute values for the *ghost nodes* are calculated depending on the values given in the tuple ```BC```. Within the tuple, we define the ```type``` (Dirichlet or Neumann) and the corresponding ```val```ue at each boundary. 

> **Note:** Usually, the results of each example within the ```GeoModBox.jl``` are stored in *gif* animations (if it is a time-dependent problem). If one wants to plot the solution for certain time steps the parameter ```save_fig``` needs to be set to 0. This setting does not result in the generation of a gif file and the single plots are not saved! Thus, care needs to be taken if the problem needs multiple time step iterations. 

> **Note:** Some examples use *named tuples* to define the different constants and variables. Alternatively, *mutable structures* can also be used to define those parameters in the ```GeoModBox.jl```. The mutable structure are especially benificial, if one needs to edit the parameters after the have been defined, e.g., for scaling. 



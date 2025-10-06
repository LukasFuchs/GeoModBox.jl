# [12 – 2D Thermal Convection (scaled)](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/12_2D_Thermal_Convection_scaled_en.ipynb)

This exercise revisits **2-D thermal convection** in a fully **nondimensional** (scaled) framework. You will define scaling constants, transform the governing equations, and study how the flow depends on the **Rayleigh number**. The setup follows the **Boussinesq approximation** with isoviscous rheology and bottom heating.

## Objectives

1. Define physically motivated **scaling constants** and apply the **nondimensional transformations** to the PDEs.  
2. Formulate the **dimensionless** energy, momentum, and mass conservation equations (incl. buoyancy term $Ra,T'$).  
3. Implement boundary conditions (Dirichlet/Neumann for temperature; **free-slip** for velocity).  
4. Solve the coupled system using finite differences (advection + diffusion + Stokes), leveraging GeoModBox solvers.  
5. Run and compare models for **$Ra = 10^4, 10^5, 10^6$**; discuss plume/slab scale and flow vigor.  
6. Compute diagnostics such as **Nusselt number** and **RMS velocity**; assess approach to steady state.

## Notes for students

- In the scaled form we assume reference parameters $\eta_0 = 1, g = 1,\kappa = 1,c_p = 1$.  
- If $Ra$ is specified, adjust $\eta_0$ accordingly; otherwise compute $Ra$ from the reference set.  
- Higher $Ra$ requires **finer grids** to maintain stability and accuracy; balance resolution vs. runtime.  
- Keep the documentation for the exercise **brief**: show the key equations, the scaling you used, and a **visualization** (animation or snapshot) of the final temperature/velocity fields.

![Exercise12a](../../assets/12_ThermalConvection_1.0e6_150_50_lineara_semilag_CNA_dc.gif)

**Figure 1.** Isoviscous, bottom-heated thermal convection for $Ra = 10^6$ with a resolution of 150x50. The initial condition is a linearly increasing temperature profile with an elliptical anomaly at the top. Thermal boundary conditions are fixed temperature at the top and bottom and zero heat flux at the sides. All velocity boundary conditions are free slip. Heat diffusion is solved using the Crank–Nicolson method, the Stokes equation using the defect correction method, and temperature advection with the semi-Lagrangian method. Models run until a steady state is reached or up to a maximum of 8000 iterations.  

![Exercise12b](../../assets/12_ThermalConvectionTimeSeries_1.0e6_150_50_lineara_semilag_CNA_dc.png)

**Figure 2.** Time series of the surface Nusselt number and the root mean square velocity. For details on how these diagnostics are calculated, see the [exercise](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/12_2D_Thermal_Convection_en.ipynb).  

![Exercise12c](../../assets/12_ThermalConvection_iterations_1.0e6_150_50_lineara_semilag_CNA_dc.png)

**Figure 3.** Variation in the root mean square velocity with numerical iterations. Empirically, a tolerance of $10^{-15}$ was chosen to define steady state. Low-$Ra$ cases typically reach steady state in fewer than 3000 iterations.  
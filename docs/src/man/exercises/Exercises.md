# Exercises

`GeoModBox.jl` contains 13 interactive exercises designed to introduce and practice the numerical solution of partial differential equations. The focus is on the **finite difference method** for solving governing equations in geodynamics.

The sequence of exercises begins with a simple 1-D Euler advection problem and proceeds to the discretization and solution of the 1-D heat diffusion equation using different numerical schemes. The heat equation is then extended to two dimensions, allowing students to analyze and compare the advantages and disadvantages of various approaches.

Next, the pure 1-D advection equation is revisited and solved with different numerical schemes, followed by coupling of the advection and diffusion equations using the operator-splitting method.

The exercises then move on to the 1-D Stokes equation, solved using both direct and iterative defect correction methods, before progressing to the steady-state and time-dependent solutions of a full Stokes problem.

The final section addresses the three fundamental conservation laws in geodynamics:  
1. Conservation of energy,  
2. Conservation of momentum, and  
3. Conservation of mass.  

These exercises investigate 2-D thermal convection, first in dimensional form and then in a scaled (non-dimensional) version.  

The series concludes with a reproduction of the **Blankenbach benchmark**.

Below, the runtimes for each exercise and their individual options are listed.  
All exercises were performed on a single CPU: *AMD Ryzen 7 7735U with Radeon Graphics*.  

| Exercise | Type | Total Runtime [s] |
|:----------|:------|:----------------|
| 1 | Euler Advection | 5.741 |
| 2 | 1-D Heat Diffusion (explicit) | 1) 1.82 s |
|   |   | 2) 1.75 s  |
|   |   | 3) 1.86 s |
| 3 | 1-D Heat Diffusion (implicit) | 1) 1.48 s  |
|   |   | 2) 1.57 s  |
|   |   | 3) 2.55 s |
| 4 | 2-D Steady State Heat Equation | 3.69 s |
| 5a | 2-D Time-dep. Heat Equation (Plume) | 1) Special Solver (Explicit): 19.74 s |
|    |   | 2) Special Solver (Implicit): 55.29 s |
|    |   | 3) General Solver (Explicit): 24.75 s |
|    |   | 4) General Solver (Implicit): 58.44 s |
|    |   | 5) General Solver (CNA): 59.35 s |
| 5b | 2-D Time-dep. Heat Equation (Sill) | 1) Special Solver (Explicit): 17.08 s |
|    |   | 2) Special Solver (Implicit): 56.61 s |
|    |   | 3) General Solver (Explicit): 31.24 |
|    |   | 4) General Solver (Implicit): 60.42 s |
|    |   | 5) General Solver (CNA): 58.9 s |
| 6 | 1-D Advection Scheme | 1) FTCS: 4.866 |
|   |   | 2) Upwind: 5.047 |
|   |   | 3) Downwind: 4.68 |
|   |   | 4) LAX: 5.202 |
|   |   | 5) SLF: 5.241 |
|   |   | 6) Semilag: 9.886 |
|   |   | 7) Tracers: 9.566 |
| 7 | 2-D Energy Equation | 1) Upwind+explicit: 9.11 s |
|   |   | 2) Upwind+implicit: 19.04 s |
|   |   | 3) Upwind+CNA: 16.35 s |
|   |   | 4) Upwind+ADI: 144.68 s |
|   |   | 5) Upwind+dc: 23.23 s |
|   |   | 6) SL+explicit: 19.1 s |
|   |   | 7) SL+implicit: 17.11 s |
|   |   | 8) SL+CNA: 16.55 s |
|   |   | 9) SL+ADI: 148.18 s |
|   |   | 10) SL+dc: 18.61 s |
|   |   | 11) Tracers+explicit: 120.41 s |
|   |   | 12) Tracers+implicit: 143.14 s |
|   |   | 13) Tracers+CNA: 143.1 s |
|   |   | 14) Tracers+ADI: 246.73 s |
|   |   | 15) Tracers+dc: 150.38 s |
| 8 | 1-D Stokes Equation | 671 ms |
| 9 | 2-D Falling Block (steady state) | 1.22 s |
| 10 | 2-D Falling Block (time-dep.) | 1) Upwind: 4.89 s |
|    |   | 2) SLF: 5.31 s |
|    |   | 3) SL: 4.94 s |
|    |   | 4) Tracers: 7.96 s |
| 11¹ | Thermal Convection (dim), **Resolution: 150×50** | **Ra = 1e4** (Diff+Adv+Momentum) |
|     |   | 1) Explicit+upwind+direct:  |
|     |   | 2) Implicit+upwind+direct:  |
|     |   | 3) CNA+upwind+direct:  |
|     |   | 4) ADI+upwind+direct:  |
|     |   | 5) DC+upwind+direct:  |
|     |   | ) Explicit+semilag+direct:  | 
|     |   | ) Implicit+semilag+direct:  | 
|     |   | ) CNA+semilag+direct:  |
|     |   | ) ADI+semilag+direct:  | 
|     |   | ) DC+semilag+direct:  | 
|     |   | ) Explicit+upwind+DC:  |
|     |   | ) Implicit+upwind+DC:  |
|     |   | ) CNA+upwind+DC:  |
|     |   | ) ADI+upwind+DC:  |
|     |   | ) DC+upwind+DC:  |
|     |   | ) Explicit+semilag+DC:  |
|     |   | ) Implicit+semilag+DC:  |
|     |   | ) CNA+semilag+DC:  |
|     |   | ) ADI+semilag+DC:  |
|     |   | ) DC+semilag+DC:  |
|     |   | **Ra = 1e5** (Diff+Adv+Momentum) |
|     |   | 1) Explicit+semilag+dc: |
|     |   | 2) CNA+semilag+dc:  |
|     |   | 3) CNA+upwind+dc:  |
|     |   | **Ra = 1e6** (Diff+Adv+Momentum)** |
|     |   | 4) Explicit+semilag+dc:  |
|     |   | 5) CNA+semilag+dc:  |
|     |   | 6) CNA+upwind+dc:  |
| 12 | Thermal Convection (scaled), **Resolution: 150×50** | **Ra = 1e4** (Diff+Adv+Momentum) |
|    |   | 1) Explicit+upwind+direct:  |
|    |   | 2) Implicit+upwind+direct:  |
|    |   | 3) CNA+upwind+direct:  |
|    |   | 4) ADI+upwind+direct:  |
|    |   | 5) DC+upwind+direct:  |
|    |   | ) Explicit+semilag+direct:  |
|    |   | ) Implicit+semilag+direct: |
|    |   | ) CNA+semilag+direct:  |
|    |   | ) ADI+semilag+direct:  |
|    |   | ) DC+semilag+direct:  |
|    |   | ) Explicit+upwind+DC:  |
|    |   | ) Implicit+upwind+DC:  |
|    |   | ) CNA+upwind+DC:   |
|    |   | ) ADI+upwind+DC: |
|    |   | ) DC+upwind+DC:   |
|    |   | ) Explicit+semilag+DC:  |
|    |   | ) Implicit+semilag+DC:  |
|    |   | ) CNA+semilag+DC:  |
|    |   | ) ADI+semilag+DC:  |
|    |   | ) DC+semilag+DC:  <br/> <br/> |
|    |   | **Ra = 1e5** (Diff+Adv+Momentum)** |
|    |   | ) Explicit+semilag+dc:  |
|    |   | ) CNA+semilag+dc:  |
|    |   | ) CNA+upwind+dc:  |
|    |   | **Ra = 1e6** (Diff+Adv+Momentum)** |
|    |   | ) Explicit+semilag+dc:  |
|    |   | ) CNA+semilag+dc:  |
|    |   | ) CNA+upwind+dc:  |
| 13 | Blankenbach Benchmark | **Resolution: 50×50** |
|    |   | Ra = 1e4:  |
|    |   | Ra = 1e5:  |
|    |   | Ra = 1e6:  |
|    |   | **Resolution: 100×100** |
|    |   | Ra = 1e4:  |
|    |   | Ra = 1e5:  |
|    |   | Ra = 1e6:  |
|    |   | **Resolution Tests:** |
|    |   | Ra = 1e4: |
|    |   | Ra = 1e5: |
|    |   | Ra = 1e6: |

¹: For the thermal convection, only the runtime where compilation and allocation overhead is avoided are shown, for the sake of simplicity. The number of iterations is shown in the parentheses. For the higher $Ra$ number models only the fastest and most accurate (see [here](../examples/Diffusion/GaussianDiffusion2D.md)) combinations of solvers from the low $Ra$ case haven been run. 
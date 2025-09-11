# Exercises

`GeoModBox.jl` contains 13 interactive exercises designed to introduce and practice the numerical solution of partial differential equations. The focus is on the **finite difference method** for solving governing equations in geodynamics.

The sequence of exercises begins with a simple 1-D Euler advection problem and proceeds to the discretization and solution of the 1-D heat diffusion equation using different numerical schemes. The heat equation is then extended to two dimensions, allowing students to analyze and compare the advantages and disadvantages of various approaches.

Next, the pure 1-D advection equation is revisited and solved with different numerical schemes, followed by coupling of the advection and diffusion equations using the operator-splitting method.

The exercises then move on to the 1-D Stokes equation, solved using both direct and iterative defect-correction methods, before progressing to the steady-state and time-dependent solutions of a full Stokes problem.

The final section addresses the three fundamental conservation laws in geodynamics:  
1. Conservation of energy,  
2. Conservation of momentum, and  
3. Conservation of mass.  

These exercises investigate 2-D thermal convection, first in dimensional form and then in a scaled (non-dimensional) version.  

The series concludes with a reproduction of the **Blankenbach benchmark**.

---

Below, the runtimes for each exercise and their individual options are listed.  
All exercises were performed on a single CPU: *AMD Ryzen 7 7735U with Radeon Graphics*.  

- The **first runtime** corresponds to the initial execution, including compilation and memory allocation.  
- The **runtime in parentheses** corresponds to a second execution, where compilation and allocation overhead are avoided.  

| Exercise Number | Type                               | Total Runtime [s] |
| --------------- | ---------------------------------- | ----------------- |
| 1               | Euler Advection                    | 5.741             |
| 2               | 1-D Heat Diffusion <br/> (explicit) | 1) 9.342 (1.14) <br/> 2) 9.975 (1.087) <br/> 3) 9.616 (1.017) |
| 3               | 1-D Heat Diffusion <br/> (implicit) | 1) 12.641 (1.749) <br/> 2) 13.537 (1.720) <br/> 3) 12.157 (1.731) | 
| 4               | 2-D Steady State Heat Equation     | 13.471 (2.069) |
| 5a              | 2-D Time-dep. Heat Equation (Plume) | 1) Explicit: 18.632 (4.799) <br/> 2) Implicit: 78.201 (60.828) |
| 5b              | 2-D Time-dep. Heat Equation (Sill)  | 1) Explicit: 17.856 (4.566) <br/> 2) Implicit: 39.626 (22.679) |
| 6               | 1-D Advection Scheme                | 1) FTCS: 11.644 (4.019) <br/> 2) Upwind: 11.03 (4.042) <br/> 3) Downwind: 11.465 (3.343) <br/> 4) LAX: 11.272 (3.405) <br/> 5) SLF: 11.312 (3.48) <br/> 6) Semilag: 14.318 (3.604) <br/> 7) Tracers: 15.662 (4.944) |
| 7               | 2-D Energy Equation                 | 1) Upwind+explicit: 17.468 (2.979) <br/> 2) Upwind+implicit: 31.441 (10.842) <br/> 3) Upwind+CNA: 32.065 (12.222) <br/> 4) Upwind+ADI: 121.665 (99.807) <br/> 5) Upwind+dc: 30.45 (11.401) <br/> 6) SLF+explicit: 22.002 (2.872) <br/> 7) SLF+implicit: 28.751 (8.911) <br/> 8) SLF+CNA: 29.876 (9.28) <br/> 9) SLF+ADI: 120.018 (98.74) <br/> 10) SLF+dc: 28.179 (9.379) <br/> 11) SL+explicit: 22.969 (6.334) <br/> 12) SL+implicit: 30.634 (9.682) <br/> 13) SL+CNA: 32.415 (11.566) <br/> 14) SL+ADI: 133.871 (110.704) <br/> 15) SL+dc: 35.364 (12.301) <br/> 16) Tracers+explicit: 89.499 (65.913) <br/> 17) Tracers+implicit: 108.464 (84.209) <br/> 18) Tracers+CNA: 112.171 (89.219) <br/> 19) Tracers+ADI: 186.71 (161.611) <br/> 20) Tracers+dc: 115.737 (89.089) |
| 8               | 1-D Stokes Equation                 | 8.489 (0.212) |



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
| 2               | 1-D Heat Diffusion\n (explicit) | 1) 9.342 (1.14) <br/> 2) 9.975 (1.087) <br/> 3) 9.616 (1.017) |

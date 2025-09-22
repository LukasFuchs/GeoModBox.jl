# [11 – 2D Thermal Convection](https://github.com/GeoSci-FFM/GeoModBox.jl/blob/main/exercises/11_2D_Thermal_Convection_en.ipynb)

This exercise introduces **2-D thermal convection** as an application of the three fundamental conservation laws:  

1. Conservation of energy,  
2. Conservation of momentum, and  
3. Conservation of mass.  

We focus on the **isoviscous, bottom-heated case**, which serves as the simplest form of mantle convection. The problem demonstrates the interplay between **advection** and **diffusion** in heat transport, while density variations due to thermal expansion drive convective flow.  

The main objectives are:  

1. Understanding the governing equations for temperature, momentum, and mass conservation,  
2. Applying the **Boussinesq approximation** to simplify the system,  
3. Introducing the **Rayleigh number** as the key nondimensional control parameter,  
4. Setting up a finite-difference discretization of the model domain,  
5. Implementing thermal and velocity boundary conditions,  
6. Solving the coupled system in time, including advection and diffusion of temperature,  
7. Analyzing convection patterns for different Rayleigh numbers ($Ra = 10^4, 10^5, 10^6$),  
8. Computing diagnostic measures such as the **Nusselt number** and RMS velocity.  

The results illustrate how increasing the Rayleigh number strengthens convection, changes the scale of plumes and slabs, and increases the overall vigor of the flow.  



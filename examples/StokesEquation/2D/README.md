# General Information

---------------------
---------------------
## Stokes Equation (2D)

## Discretization 

&emsp;The conservation equations of *momentum* and *mass* are solved properly in two dimensions (*x* and *y*) using a staggered finite difference grid, where the horizontal and vertical velocity are defined in between the regular grid points, and the pressure within a finite difference cell (Figure 1). A staggered grid enables the conservation of the stress between adjacent grid points and one can solve equations (8) and (9) for the unknows.  

<img src="./Figures/MomentumGrid_2D.png" alt="drawing" width="600"/> <br>
**Figure 1. Staggered finite difference grid.** Discretization of the conservation equations of momemtum and mass. The horizontal and vertical velocities are defined in between the vertices (cyan and orange lines, respectively), and the pressure is defined on the centroids. The horizontal and vertical velocities require *ghost nodes* at the north and south and east and west boundary, respectively. 

### Solution 

### Boundary Conditions

-----------------------
-----------------------

## Examples

<!-- 
- 1D case 
-- Discretized equations
-- Solving the equations 
--- Direct solution 
--- Defection corrections solution
- 2D case 
-->

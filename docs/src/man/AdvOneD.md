# Advection Equation (1D)

$\begin{equation}
\frac{\partial{T}}{\partial{t}} = -v_x \left(\frac{\partial{T}}{\partial{x}}\right),
\end{equation}$

- Advection equation in 1D 

- *FTCS* alway instable, can be tested within the exercises

$\begin{equation}
T_i^{n+1} = T_i^n - v_x \Delta{t}\frac{T_{i+1}^n - T_{i-1}^n}{2\Delta{x}},
\end{equation}$

- Upwind -> Courtant criterium, numerische diffusion, exact solution for \alpha = 1
    -> error estimaiton of each FD term, taylor expansion 
    -> \alpha = 1 no numerical diffusion, resolution dependent

$\begin{equation}
\frac{T_{i}^{n+1}-T_{i}^n}{\Delta{t}} = -v_{x,i}
\end{equation}$

- Lax Method -> stable but strong numerical diffusion

$\begin{equation}
\end{equation}$

- Staggered Leap Frog -> Higher Fehlerordnung in time and space, no numerical diffusion 

$\begin{equation}
\end{equation}$

- Semi-lagrangian

$\begin{equation}
\end{equation}$

- Tracer

$\begin{equation}
\end{equation}$


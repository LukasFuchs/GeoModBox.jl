# [Oceanic Geotherm](https://github.com/LukasFuchs/GeoModBox.jl/blob/main/examples/DiffusionEquation/1D/OceanicGeotherm_1D.jl) 

The 1-D temperature profile of an oceanic geotherm can be calculated by solving the conductive part of the 1-D *temperature equation* using variable thermal parameters with a proper conserving finite difference scheme (so far only including a radiogenic heat source). For the sake of continuity, we use the 1-D solver for variable thermal parameters, even though, we assume a constant thermal conductivity in this example.

A proper conservative finite difference scheme in 1-D means that the temperature is defined on the *centroids* and the vertical heat flux and the thermal conductivity $k$ on the *vertices*.

The 1-D temperature equation is given by: 

$\begin{equation}
\rho c_{p} \frac{\partial{T}}{\partial{t}} = \frac{\partial{q_y}}{\partial{y}} + \rho H,
\end{equation}$ 

where the heat flux $q_y$ is defined as:

$\begin{equation}
\left. q_{y,m} = -k_m \frac{\partial T}{\partial y}\right\vert_{m},\ \textrm{for}\ m = 1:nv, 
\end{equation}$

with $\rho, c_{p}, T, t, k, H, y$ and $nv$ are the density [kg/m³], the specific heat capacity [J/kg/K], the temperature [K], the time [s], the thermal conductivity [W/m/K], the heat generation rate per mass [W/kg], the depth [m], and the number of vertices, respectively. 

For more details on how to discretize the equation using an explicit, forward Euler finite difference scheme see the [documentation](./DiffOneD.md).

> **Note:** Here, we use *named tuples* to define the different constants and variables. Alternatively, mutable structures to define those parameters are also available in the ```GeoModBox.jl```. 

Let's start with the definition of the geometrical, numerical, and physical constants: 

```Julia 
# Constants --------------------------------------------------------- #
H           =   200e3               #   Hight of the model [ m ]
nc          =   200                 #   Number of central grid points    
Δy          =   H/nc                #   Grid resolution

# Depth ---
yc          =   LinRange(-H + Δy/2.0,0.0 - Δy/2.0,nc)     
yv          =   LinRange(-H,0,nc+1)
    
Py  =   (
    ρm      =   3000,               #   Density [ kg/m^3 ]
    cpm     =   1000,               #   Heat capacity [ J/kg/K ]
    km      =   3.0,                #   Conductivity [ W/m/K ]
    HM      =   0,                  #   Heat generation rate [W/kg]; Q = ρ*H0
)    
# ------------------------------------------------------------------- #
```

In the following, we need to define the initial and boundary condition: 

1. Temperature at the surface and bottom
2. Linear increasing temperature profile assuming a certain adiabatic gradient and potential mantle temperature


```Julia
# Initial Condition ------------------------------------------------- #
T   =   (
    Tpot    =   1315 + 273.15,      #   Potential temperautre [ K ]
    ΔTadi   =   0.5,                #   Adiabatic temperature gradient [ K/km ]
    Ttop    =   273.15,             #   Surface temperature [ K ]
    T_ex    =   zeros(nc+2,1),    
)
T1  =   (
    Tbot    =   T.Tpot + T.ΔTadi*abs(H/1e3),    # Bottom temperature [ K ]
    T       =   T.Tpot .+ abs.(yc./1e3).*T.ΔTadi,   # Initial T-profile [ K ]
)
T   =   merge(T,T1)
     
Tini                =   zeros(nc,1)
Tini                .=   T.T
T.T_ex[2:end-1]     .=   T.T
# ------------------------------------------------------------------- #
```

We can choose between *Dirichlet* or *Neumann* thermal boundary conditions at the surface and the bottom. The boundary conditions are explicitly implemented in the solver, where the temperature for the *ghost nodes* is calculated depending on the values given in the tuple ```BC```. Within the tuple, we define the ```type``` and the corresponding ```val```ue of each boundary (North or South). 

```Julia 
# Boundary conditions ----------------------------------------------- #
BC      =   (
    type    = (N=:Dirichlet, S=:Dirichlet),
    val     = (N=T.Ttop[1],S=T.Tbot[1])
)
# If Neumann boundary conditions are choosen, the following values result in the given heatflux for the given thermal conductivity k. 
# S      =   -0.03;          # c     =   -k/q -> 90 mW/m^2
# N      =   -0.0033;        # c     =   -k/q -> 10 mW/m^2
# ------------------------------------------------------------------- #
```

In the following, we need to define the multiplication factor ```fac``` of the *diffusion stability criterion* for the explicit thermal solver. 

```Julia
# Time stability criterion ------------------------------------------ #
fac     =   0.8                 #   Courant criterion
tmax    =   60                  #   Lithosphere age [ Ma ]
tsca    =   60*60*24*365.25     #   Seconds per year

age     =   tmax*1e6*tsca        #   Age in seconds    
# ------------------------------------------------------------------- #
```

Let's check that the initial condition is properly defined. 

```Julia
# Plot Initial condition -------------------------------------------- #
plotparam   =   1
q = plot(Tini,yc./1e3, 
    label="", 
    xlabel="T [ K ]", ylabel="z [ km ]", 
    title="Initial Temperature",
    xlim=(T.Ttop,T.Tbot),ylim=(-H/1e3,0))
display(q)
# ------------------------------------------------------------------- #
```

![OG1D_iniT](../assets/OG1D_iniT.svg)

**Figure 1. Initial temperature profile.**

Since we use the thermal solver for variable thermal parameters, we need to expand the scalar to a vector with the dimensions of the number of centroids ```nc```. Additionally, we define the thermal diffusivity $\kappa$ and initialize the vertical heat flux ```q```. 

```Julia
# Setup Fields ------------------------------------------------------ #
Py1     =   (
    ρ       =   Py.ρm.*ones(nc,1),
    cp      =   Py.cpm.*ones(nc,1),
    k       =   Py.km.*ones(nc+1,1),
    H       =   Py.HM.*ones(nc,1)
)
Py  =   merge(Py,Py1)
Py2     =   (
    # Thermal diffusivity [ m^2/s ] 
    κ       =  maximum(Py.k)/minimum(Py.ρ)/minimum(Py.cp),     
)
Py  =   merge(Py,Py2)
T2  =   (
    q   =   zeros(nc+1,1),
)
T   =   merge(T,T2)
# ------------------------------------------------------------------- #
```

...


```Julia
# Time stability criterion ------------------------------------------ #
Δtexp   =   Δy^2/2/Py.κ             #   Stability criterion for explicit
Δt      =   fac*Δtexp               #   total time step

nit     =   ceil(Int64,age/Δt)      #   Number of iterations    

time    =   zeros(1,nit)            #   Time array
# ------------------------------------------------------------------- #
```

```Julia
# Calculate 1-D temperature profile --------------------------------- #
count   =   1
for i = 1:nit
    if i > 1
        time[i]     =   time[i-1] + Δt
    elseif time[i] > age
        Δt          =   age - time[i-1]
        time[i]     =   time[i-1] + Δt
    end
    ForwardEuler1D!(T,Py,Δt,Δy,nc,BC)
    if i == nit || abs(time[i]/1e6/tsca - count*5.0) < Δt/1e6/tsca        
        println(string("i = ",i,", time = ", time[i]/1e6/tsca))        
        q = plot!(T.T,yc./1e3, 
            label=string("t = ",ceil(time[i]/1e6/tsca),"[Ma]"), 
            xlabel="T [ K ]", ylabel="z [ km ]", 
            title="Oceanic Geotherm",
            xlim=(T.Ttop,T.Tbot),ylim=(-H/1e3,0))
        display(q)
        count = count + 1
    end
end
# ------------------------------------------------------------------- #
```

... 

```Julia
# Calculate heaf flow ----------------------------------------------- #
# South ---
#T.T_ex[2:end-1]     =   T.T
T.T_ex[1]   =   (BC.type.S==:Dirichlet) * (2 * BC.val.S - T.T_ex[2]) + 
                (BC.type.S==:Neumann) * (T.T_ex[2] - BC.val.S*Δy)
# North ---
T.T_ex[end] =   (BC.type.N==:Dirichlet) * (2 * BC.val.N - T.T_ex[nc+1]) +
                (BC.type.N==:Neumann) * (T.T_ex[nc+1] + BC.val.N*Δy)
if size(Py.k,1)==1
    for j=1:nc+1
        T.q[j]  =   -Py.k * 
            (T.T_ex[j+1] - T.T_ex[j])/Δy
    end    
else
    for j=1:nc+1
        T.q[j]  =   -Py.k[j] * 
            (T.T_ex[j+1] - T.T_ex[j])/Δy
    end
end
# ------------------------------------------------------------------- #
```

...

```Julia
# Plot -------------------------------------------------------------- #
if BC.type.N==:Dirichlet && BC.type.S==:Dirichlet
    Tana    =   zeros(nc,1)
    Tana    .=   Tini .+ (T.Ttop - T.Tpot).*erfc.(-yc./(2*sqrt(age*Py.κ)))
    Tana[end] =   T.Ttop
end    
p = plot(T.T,yc./1e3, 
        label=string("t = ",ceil(maximum(time)/1e6/tsca),"[Ma]"), 
        xlabel="T [ K ]", ylabel="z [ km ]",
        title="Oceanic Geotherm",
        xlim=(T.Ttop,T.Tbot),ylim=(-H/1e3,0),
        layout=(1,2),subplot=1)        
if BC.type.N==:Dirichlet && BC.type.S==:Dirichlet
    plot!(p,Tana,yc./1e3, 
            label="T_HSCM",linestyle=:dash,
            layout=(1,2),subplot=1)        
end        
p = plot!(T.q.*1e3,yv./1e3, 
        label="", 
        xlabel="q [ mW ]", ylabel="z [m]", 
        title="Heat Flux",
        ylim=(-H/1e3,0),
        subplot=2)        
display(p)
savefig(p,"./examples/DiffusionEquation/1D/Results/OceanicGeotherm_1D.png")
savefig(q,"./examples/DiffusionEquation/1D/Results/OceanicGeotherm_1D_evolve.png")
# ======================================================================= #
```

...

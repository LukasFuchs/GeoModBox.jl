using Plots, SpecialFunctions
using GeoModBox.HeatEquation.OneD

@doc raw"""
    OceanicGeotherm_1D()

Function to calculate the 1D geotherm for an oceanic lithosphere. 
Temperature is calculated by solving the 1-D heat equation assuming
variable thermal parameters and a radiogenic heat source.
The equation is solved using a proper conserving finite difference
scheme.

"""
function OceanicGeotherm_1D()
# ------------------------------------------------------------------- #
#    LF - 22.10.2024 - juila                                          #
# ------------------------------------------------------------------- #
    
# Constants --------------------------------------------------------- #
H           =   200e3               #   Hight of the model [ m ]
nc          =   200                 #   Number of central grid points    
Δy          =   H/nc                #   Grid resolution

# Depth [ m ] ---
yc          =   LinRange(-H + Δy/2.0,0.0 - Δy/2.0,nc)     
yv          =   LinRange(-H,0,nc+1)
    
Py  =   (
    ρm      =   3000,               #   Density [ kg/m^3 ]
    cpm     =   1000,               #   Heat capacity [ J/kg/K ]
    km      =   3.0,                #   Conductivity [ W/m/K ]
    HM      =   0,                  #   Heat generation rate [W/kg]; Q = ρ*H0
)    
# ------------------------------------------------------------------- #
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
# Boundary conditions ----------------------------------------------- #
BC      =   (
    type    = (N=:Dirichlet, S=:Dirichlet),
    val     = (N=T.Ttop[1],S=T.Tbot[1])
)
# S      =   -0.03;          # c     =   -k/q -> 90 mW/m^2
# N      =   -0.0033;        # c     =   -k/q -> 10 mW/m^2
# ------------------------------------------------------------------- #
# Time stability criterion ------------------------------------------ #
fac     =   0.8                 #   Courant criterion
tmax    =   60                  #   Lithosphere age [ Ma ]
tsca    =   60*60*24*365.25     #   Seconds per year

age     =   tmax*1e6*tsca        #   Age in seconds    
# ------------------------------------------------------------------- #
# Plot Initial condition -------------------------------------------- #
plotparam   =   1
q = plot(Tini,yc./1e3, 
    label="", 
    xlabel="T [ K ]", ylabel="z [ km ]", 
    title="Initial Temperature",
    xlim=(T.Ttop,T.Tbot),ylim=(-H/1e3,0))
display(q)
# ------------------------------------------------------------------- #
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
# Time stability criterion ------------------------------------------ #
Δtexp   =   Δy^2/2/Py.κ             #   Stability criterion for explicit
Δt      =   fac*Δtexp               #   total time step

nit     =   ceil(Int64,age/Δt)      #   Number of iterations    

time    =   zeros(1,nit)            #   Time array
# ------------------------------------------------------------------- #
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
end

OceanicGeotherm_1D()
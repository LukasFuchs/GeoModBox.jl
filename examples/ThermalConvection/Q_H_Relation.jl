using Printf

function Q_H_Relation()

# Scaling parameters ---------------------------------------------------- #
k       =   4.125           # [ W/m/K ]
ρ       =   3300            # [ kg/m^3 ]
cp      =   1250            # [ J/kg/K ]
ΔT      =   2500            # [ K ]

κ       =   1e-6            # [ m^2/s ]
g       =   9.81            # [ m/s^2 ] 
α       =   2.0e-5          # [ K^-1 ]

# R       =   6.371e6         # [ m ]
D       =   2.900e6         # [ m ] 2.871e6
ar      =   4               # Aspect ratio
L       =   ar * D          # Length

# Rcmb    =   R-D             # [ m ]

η       =   3.947725485e22  # [ Pa s ]

RaD     =   ρ*g*α*ΔT*D^3/(η*κ)

@printf("- Calculate Heat Production Rate -")
@printf("\n")
@printf(" RaD: %g\n",RaD)

# ## Mantle Volume: V[M] = V[R] - V[C] = 4/3 * pi * [ R^3 - (R-D)^3 ]
# #V_M     =   4/3 *π * ( R^3 - (R-D)^3 )
# # Mantle Area: L*D [ m² ]
# V_A     =   L*D

# # @printf(" Mantle Volume [m^3]: %g\n ",V_M)
# @printf(" Mantle Volume [m^2]: %g\n ",V_A)

# # Earth Surface Area (O = 4 * pi * R^2): 
# O       =   4 * pi * R^2

# # Average surface heat flux ( 90 mW/m^2 ) ------------------------------- #
# @printf("\n-----------------------\n")
# @printf(" EARTH DATA \n")
# @printf("")
# q_mean  = 90e-3

# # Total energy loss: 
# E_total =   q_mean * O
# @printf(" Total Energyloss (Earth - 90 mW/m^2) [TW]: %2.3g\n",
#                 E_total/1e12)

# Heat generation rate -------------------------------------------------- #
@printf("\n-----------------------\n")
@printf(" MODEL DATA \n")
Qp      =   25; 

Qsc     =   (ΔT*κ*ρ*cp) / D^2
Hsc     =   (ΔT * k)/ (D^2 * ρ)

# Heat generation per mass (dimensional): 
Q       =   Qp * Qsc 
H       =   Q / ρ
Hp      =   H / Hsc

@printf(" Heat Production Rate per Volume (non-dim): %g\n",Qp)
@printf(" Heat Production Rate per mass (non-dim): %g\n\n",Hp)

@printf(" Heat generation rate per Volume [W/m^3]: %1.2e\n",Q)
@printf(" Heat generation rate per mass [W/kg]: %1.2e\n",H)


# Q_total =   Q * V_M
# @printf(" Total Heatenergie (Model - H = %g) [TW]: %2.3g\n",
#                 Hp,Q_total/1e12)

# q_model =   Q_total / O
# @printf(" Average Surface Heat Flux (Model - H = %g) [mW/m^2]: %2.3g\n",
#                 Hp,q_model*1e3)

# qsc     =   k * DeltaT / R; 
            
# q_modsc =   q_model / qsc;
# @printf([" Average Surface Heat Flux (Model - H = ",num2str(Hp),", non-dim): ",...
#                 sprintf("#2.3g\n",q_modsc)])
            
# H_mod   =   q_mean * O / (rho * V_M);
# H_modp  =   H_mod / Hsc; 
# @printf([" Modified Heat Production Rate (non-dim): ",num2str(H_modp)])
# @printf(" (to get mean Earth surface heat flux for purely internally heated model)")

# Ra_range = [1e5,1e6,1e7,1e8,1e9];
# T_i     =   0.5; 
# H_ideal =   (T_i.*Ra_range.^(1/4)).^(4/3)


# # Scaling relationship for mixed heating: T ~ Q^0.75/Ra^0.25;
# Q1  = Hp; 
# Ra1 = RaD;

# Ra2 = 1e7; 
# Q2  = Q1*(Ra2/Ra1)^(1/3)

end

Q_H_Relation()
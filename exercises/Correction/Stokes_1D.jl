using Plots, ExtendableSparse

function Stokes_1D_const()
# Model Parameter ------------------------------------------------------- #
M   =   (
    ymin        =   -400.0e3,           #   Tiefe [ m ]
    ymax        =   0.0e3,              
)
I   =   (
    vₓ₀         =   5/100/31536000,     #   Geschwindigkeit oben [ m/s ]
    η₀          =   1.0e21,             #   Viscosität [ Pa s ]
    ∂P∂x        =   -5.0e1,              #   horizontaler Druckgradient [ Pa/m ]
)
# ----------------------------------------------------------------------- #
# Numerische Parameter -------------------------------------------------- #
NC  =   (
    y   =   100,        #   Anzahl der vertikalen Centroids
)
Δ   =   (
    y   =   (M.ymax-M.ymin)/NC.y,   #   Gitterauflösung
)
# ----------------------------------------------------------------------- #
# Speicherzuweisung ----------------------------------------------------- #
D   =   (
    vₓ  =   zeros(NC...),
    vₓₐ =   zeros(NC...),
    Δvₓ =   zeros(NC...)
)
# ----------------------------------------------------------------------- #
# Gitter ---------------------------------------------------------------- #
y   =   (
    c   =   LinRange(M.ymin+Δ.y/2,M.ymax-Δ.y/2,NC.y),
)
# ----------------------------------------------------------------------- #
# Analytische Lösung ---------------------------------------------------- # 
@. D.vₓₐ    =   1.0/2.0/I.η₀ * I.∂P∂x * 
                    (y.c^2 + (M.ymax-M.ymin).*y.c) + 
                    I.vₓ₀*y.c/(M.ymax-M.ymin) + I.vₓ₀
# ----------------------------------------------------------------------- #
# Randbedingungen ------------------------------------------------------- #
VBC     =   (
    type    = (S=:Dirichlet, N=:Dirichlet),
    val     = (S=0.0,N=I.vₓ₀)
)
# ----------------------------------------------------------------------- #
# Lineares Gleichungssystem --------------------------------------------- #
Num     =   (Vx=1:NC.y,)
ndof    =   maximum(Num.Vx)
K       =   ExtendableSparseMatrix(ndof,ndof)
rhs     =   zeros(NC...)
# ----------------------------------------------------------------------- #
# Zusammenstellen der Koeffizientenmatrix ------------------------------- #
a       =     I.η₀ / Δ.y^2.0
b       =   - 2.0 * I.η₀ / Δ.y^2.0
rhs     .=  I.∂P∂x
for i = 1:NC.y
    # Gleichungsnummer ---
    ii  =   Num.Vx[i]
    # Stempel ---
    iN  =   ii + 1      #   Norden
    iC  =   ii          #   Zentral
    iS  =   ii - 1      #   Süden
    # Ränder ---
    # Falls ein Süd Index gebrauch wird ---
    inS    =  i==1    ? false  : true
    DirS   = (i==1    && VBC.type.S==:Dirichlet) ? 1. : 0.
    NeuS   = (i==1    && VBC.type.S==:Neumann  ) ? 1. : 0.
    # If an East index is required ---
    inN    =  i==NC.y ? false  : true
    DirN   = (i==NC.y && VBC.type.N==:Dirichlet) ? 1. : 0.
    NeuN   = (i==NC.y && VBC.type.N==:Neumann  ) ? 1. : 0.
    if inS K[ii,iS]    = a end
        K[ii,iC]       =   b + (NeuS + NeuN - DirS - DirN)*a
    if inN K[ii,iN]    = a end    
    # Änderung der rechten Seite ---
    rhs[i]      +=  a*VBC.val.S*Δ.y * NeuS - 
                        2*a*VBC.val.S * DirS - 
                        a*VBC.val.N*Δ.y * NeuN - 
                        2*a*VBC.val.N * DirN
end
# ----------------------------------------------------------------------- #
# Lösung des Gleichungssystems ------------------------------------------ #
D.vₓ      .=   K \ rhs
# ----------------------------------------------------------------------- #
# Abweichung vom der analytischen Lösung -------------------------------- #
@. D.Δvₓ    =   ((D.vₓₐ - D.vₓ) / D.vₓₐ) * 100.0
# ----------------------------------------------------------------------- #
# Visualisierung -------------------------------------------------------- #
q   =   plot(D.vₓ,y.c./1e3,label="numerical",
            xlabel="vₓ", ylabel="y [km]",
            title="Velocity Profile",
            xlim=(0,I.vₓ₀*1.5),ylim=(M.ymin/1e3,M.ymax/1e3),
            layout=(1,2),subplot=1)
plot!(q,D.vₓₐ,y.c./1e3,label="analytical",linestyle=:dash,
        subplot=1)
plot!(q,D.Δvₓ,y.c./1e3,label="",
        xlabel="Δvₓ",ylabel="z [km]",
        title="Error [ % ]",
        subplot=2)
display(q)
# ----------------------------------------------------------------------- #
end

Stokes_1D_const()











# clear
# clc

# ## Model Geometry ======================================================= #
# H       =   -400e3;             # [ m ]
# vx0     =   2;                  # [ cm/a ]
# vx0     =   vx0/100/31536000;   # [ m/s ]
# dPdx    =   -1e2;               # [ Pa/m ]; z.B. -7e1 (Warum negativ?)

# ## Numerische Parameter ================================================= #
# nz      =   51;
# dz      =   H/(nz-1);

# # Koordinaten
# z       =   H:-dz:0;
# # Index der internen Gitterpunkte
# indi    =   2:nz-1;

# ## Definition der Vikosität ============================================= #
# eta0        =   22;                 # Power der Viskosität oben
# eta1        =   21;                 # Power der Viskosität unten
# m           =   10^eta1/10^eta0;    # Viskositätsverhältnis
# # eta         =   10^eta0.*m.^(z'./H);   # Viskositätsprofil
# eta         =   10.^eta0*exp(log(m).*z'./H);

# ## Analytische Loesung ================================================== #
# vx_ana      =   -dPdx*H/10^eta0/log(m)/(m-1).*...
#     (z.*(m.^((-z+H)./H) - m.^(-z./H)) + H.*(m.^(-z./H) - 1)) + ...
#     vx0/(m-1).*(m.^((-z+H)./H) - 1);

# ## Erstellung der Koeffizientenmatrix A
# diag        =   zeros(nz,3);

# # Bestimmung der Diagonalen
# diag(indi-1,1)  =   (eta(indi)+eta(indi-1))/2/dz^2;
# diag(indi,2)    =   -(eta(indi-1)+2*eta(indi)+eta(indi+1))/2/dz^2;
# diag(indi+1,3)  =   (eta(indi)+eta(indi+1))/2/dz^2;

# # Randbedingungen - no slip, d.h. konstante Geschwindigkeit
# diag(1,2)       =   1;
# diag(nz,2)      =   1;

# A               =   spdiags(diag,[-1 0 1],nz,nz);

# ## Definition der rechten Seite des Gleichungssystems
# rhs             =   zeros(nz,1);

# rhs(indi)       =   dPdx;
# rhs(1)          =   0;          # unten
# rhs(nz)         =   vx0;        # oben

# vx              =   A\rhs;

# ## Darstellung der Daten
# figure(1)
# clf
# subplot(1,2,1)
# plot(vx,z./1e3,'LineWidth',2)
# hold on
# plot(vx_ana,z./1e3,'r--','LineWidth',2)
# xlabel('v_x [ m/s ]'); ylabel('z [km]'); title('Geschwindigkeitsprofil')
# set(gca,'FontWeight','Bold','FontSize',15,'LineWidth',2);
# subplot(1,2,2)
# plot(eta,z./1e3,'LineWidth',2)
# xlabel('\eta[ Pa s ]'); ylabel('z [km]'); title('Viskosity')
# set(gca,'xscale','log','FontWeight','Bold','FontSize',15,'LineWidth',2);


# ## ===================================================================== ##
# # ============================== ENDE =================================== #
# # ======================================================================= #













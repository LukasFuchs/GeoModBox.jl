using Plots, ExtendableSparse

function Stokes_1D_const()
# Model Parameter ------------------------------------------------------- #
M   =   (
    ymin        =   -400.0e3,           #   Tiefe [ m ]
    ymax        =   0.0e3,              
)
I   =   (
    vₓ₀         =   5/100/31536000,     #   Geschwindigkeit oben [ m/s ]
    η₀          =   1.0e21,             #   Viscosität an der Oberseite [ Pa s ]
    η₁          =   1.0e21,             #   Viskosität an der Unterseite [ Pa s ]
    ∂P∂x        =   -5.0e1,              #   horizontaler Druckgradient [ Pa/m ]
)
I1  =   (
    m           =   I.η₁ / I.η₀,        #   Viskositätsverhältnis
)
I   =   merge(I,I1)
# ----------------------------------------------------------------------- #
# Numerische Parameter -------------------------------------------------- #
NC  =   (
    y   =   100,        #   Anzahl der vertikalen Centroids
)
NV  =   (
    y   =   NC.y + 1,
)
Δ   =   (
    y   =   (M.ymax-M.ymin)/NC.y,   #   Gitterauflösung
)
# ----------------------------------------------------------------------- #
# Speicherzuweisung ----------------------------------------------------- #
D   =   (
    η   =   zeros(NC.y+1),
    vₓ  =   zeros(NC...),
    vₓₐ =   zeros(NC...),
    Δvₓ =   zeros(NC...)
)
# ----------------------------------------------------------------------- #
# Gitter ---------------------------------------------------------------- #
y   =   (
    c   =   LinRange(M.ymin+Δ.y/2,M.ymax-Δ.y/2,NC.y),
    v   =   LinRange(M.ymin,M.ymax,NV.y),
)
@. D.η  =   I.η₀ * exp(-log(I.m)* y.v / (M.ymax-M.ymin))
# ----------------------------------------------------------------------- #
# Analytische Lösung ---------------------------------------------------- #
if I.m  == 1.0
    @. D.vₓₐ    =   1.0/2.0/I.η₀ * I.∂P∂x * 
                    (y.c^2 + (M.ymax-M.ymin).*y.c) + 
                    I.vₓ₀*y.c/(M.ymax-M.ymin) + I.vₓ₀
else
    @. D.vₓₐ    =   -I.∂P∂x * (M.ymax-M.ymin) / I.η₀ / log(I.m) / (I.m-1) * 
        (-y.c * (I.m^((y.c + (M.ymax-M.ymin))/(M.ymax-M.ymin)) - 
        I.m^(y.c/(M.ymax-M.ymin))) + (M.ymax-M.ymin)*(I.m^(y.c/(M.ymax-M.ymin)) - 1)) + 
        I.vₓ₀ / (I.m-1) * (I.m ^ ((y.c+(M.ymax-M.ymin))/(M.ymax-M.ymin)) - 1)
end
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
#a       =     I.η₀ / Δ.y^2.0
#b       =   - 2.0 * I.η₀ / Δ.y^2.0
rhs     .=  I.∂P∂x
for i = 1:NC.y
    a   =   D.η[i] / Δ.y^2.0
    b   =   -(D.η[i]+D.η[i+1]) / Δ.y^2.0
    c   =   D.η[i+1] / Δ.y^2.0
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
        K[ii,iC]       =   b + (NeuS - DirS)*a + (NeuN - DirN)*c
    if inN K[ii,iN]    = c end    
    # Änderung der rechten Seite ---
    rhs[i]      +=  a*VBC.val.S*Δ.y * NeuS - 
                        2*a*VBC.val.S * DirS - 
                        c*VBC.val.N*Δ.y * NeuN - 
                        2*c*VBC.val.N * DirN
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
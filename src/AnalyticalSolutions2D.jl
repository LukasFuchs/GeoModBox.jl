function Dani_Solution_vec!(type,D,M,x,y,rad,mus_i,NC,NV)
# [ Vx_N,Vx_S,Vx_W,Vx_E,Vy_N,Vy_S,Vy_W,Vy_E,Pa,Vxa,Vya ]
# ------------------------------------------------------------------------ #
# ANALYTICAL SOLUTION - PRESSURE AND VELOCITY AROUND A CIRCULAR INCLUSION
# BASED ON DANI SCHMID'S 2002 CYL_P_MATRIX.M
# Vectorised version by:
# Thibault Duretz, Ludovic Raess - Unil 2016
# ------------------------------------------------------------------------ #
    if type==:SimpleShear
        gr  =   1.0
        er  =   0.0
    elseif type==:PureShear
        gr  =   0.0
        er  =   -1.0
    end
    mm      =   1.0               #     Viscosity of matrix
    mc      =   mus_i
    A       =   mm.*(mc-mm)./(mc+mm)
    i       =   im
    # -------------------------
    # Pa      =   zeros(Float64,NC...)
    # Vxa     =   zeros(Float64,NV.x,NC.y)
    # Vya     =   zeros(Float64,NC.x,NV.y)

    # PRESSURE ----------------------------------------------------------- #
    XC2     =   copy(x.c2d) .- (M.xmax-M.xmin)/2 # x.c .+ 0*y.c'
    YC2     =   copy(y.c2d) .+ (M.ymax-M.ymin)/2 # 0*x.c .+ y.c'
    Z       =   zeros(NC...)
    Z       =   XC2 .+ i .* YC2
    PH      =   falses(NC...)
    PH      .=  (XC2.^2 .+ YC2.^2 .<= rad^2)
    P       =   -2 * mm * (mc - mm) / (mc + mm) * real.(rad^2 ./ Z.^2 .* (i * gr + 2 * er))
    D.Pa[.!PH]  .=  P[.!PH]

    # Conforming Nodes --------------------------------------------------- #
    # VELOCITY X ---
    XV2     =   copy(x.vx2d[:,2:end-1]) .- (M.xmax-M.xmin)/2
    YC2     =   copy(y.vx2d[:,2:end-1]) .+ (M.ymax-M.ymin)/2
    Z       =   zeros(NV.x,NC.y)
    Z       =   XV2 .+ i .* YC2
    PH      =   falses(NV.x, NC.y)
    PH      .=  (XV2.^2 .+ YC2.^2 .<= rad^2)
    V_tot   =   (mm / (mc + mm)) * (i * gr + 2 * er) .* conj.(Z) .- (i / 2) * gr .* Z
    D.Vxa[PH] .=  real.(V_tot[PH])

    phi_z   =   -(i/2) * mm * gr .* Z .- (i * gr + 2 * er) * A * rad^2 ./ Z
    d_phi_z =   -(i/2) * mm * gr .+ (i * gr + 2 * er) * A * rad^2 ./ Z.^2
    psi_z   =   (i * gr - 2 * er) * mm .* Z .- (i * gr + 2 * er) * A * rad^4 ./ Z.^3

    V_tot   =   (phi_z .- Z .* conj.(d_phi_z) .- conj.(psi_z)) ./ (2 * mm)
    D.Vxa[.!PH]   .=  real.(V_tot[.!PH])

    # VELOCITY Y ---
    XC2     =   copy(x.vy2d[2:end-1,:]) .- (M.xmax-M.xmin)/2
    YV2     =   copy(y.vy2d[2:end-1,:]) .+ (M.ymax-M.ymin)/2
    Z       =   zeros(NC.x,NV.y)
    Z       =   XC2 .+ i .* YV2
    PH      =   falses(NC.x,NV.y)
    PH      .=  (XC2.^2 .+ YV2.^2 .<= rad^2)
    V_tot   =   (mm / (mc + mm)) * (i * gr + 2 * er) .* conj.(Z) .- (i / 2) * gr .* Z
    D.Vya[PH] .=  imag.(V_tot[PH])

    phi_z   =   -(i/2) * mm * gr .* Z .- (i * gr + 2 * er) * A * rad^2 ./ Z
    d_phi_z =   -(i/2) * mm * gr .+ (i * gr + 2 * er) * A * rad^2 ./ Z.^2
    psi_z   =   (i * gr - 2 * er) * mm .* Z .- (i * gr + 2 * er) * A * rad^4 ./ Z.^3

    V_tot   =   (phi_z .- Z .* conj.(d_phi_z) .- conj.(psi_z)) ./ (2 * mm)
    D.Vya[.!PH] .= imag.(V_tot[.!PH])

    # Get BC ---
    D.Vx_W[:]    =   D.Vxa[1   , :  ]'
    D.Vx_E[:]    =   D.Vxa[end , :  ]'
    D.Vy_S[:]    =   D.Vya[:   ,1   ]
    D.Vy_N[:]    =   D.Vya[:   ,end ]
    
    # Non Conforming Nodes ----------------------------------------------- #
    Vx_NC   =   zeros(Float64,NV...)
    Vy_NC   =   zeros(Float64,NV...)
    # Vx_S    =   zeros(Float64,NV.x,1)
    # Vx_N    =   zeros(Float64,NV.x,1)
    # Vy_W    =   zeros(Float64,1,NV.y)'
    # Vy_E    =   zeros(Float64,1,NV.y)'

    # VELOCITY X & Y ---
    XV2     =   copy(x.v2d) .- (M.xmax-M.xmin)/2
    YV2     =   copy(y.v2d) .+ (M.ymax-M.ymin)/2
    Z       =   copy(YV2)
    Z       =   XV2 .+ i .* YV2
    PH      =   falses(NV.x,NV.y)
    PH      .=  (XV2.^2 .+ YV2.^2 .<= rad^2)
    V_tot   =   (mm / (mc + mm)) * (i * gr + 2 * er) .* conj.(Z) .- (i / 2) * gr .* Z
    Vx_NC[PH]   .=  real.(V_tot[PH])
    Vy_NC[PH]   .=  imag.(V_tot[PH])

    phi_z   =   -(i/2) * mm * gr .* Z .- (i * gr + 2 * er) * A * rad^2 ./ Z
    d_phi_z =   -(i/2) * mm * gr .+ (i * gr + 2 * er) * A * rad^2 ./ Z.^2
    psi_z   =   (i * gr - 2 * er) * mm .* Z .- (i * gr + 2 * er) * A * rad^4 ./ Z.^3

    V_tot   =   (phi_z .- Z .* conj.(d_phi_z) .- conj.(psi_z)) ./ (2 * mm)
    Vx_NC[.!PH]     .=  real.(V_tot[.!PH])
    Vy_NC[.!PH]     .=  imag.(V_tot[.!PH])

    # Get BC ---
    D.Vx_S[:,1] = Vx_NC[:,1  ]
    D.Vx_N[:,1] = Vx_NC[:,end]
    D.Vy_W[:,1] = Vy_NC[1  ,:]'
    D.Vy_E[:,1] = Vy_NC[end,:]'
    
    # return Vx_N, Vx_S, Vx_W, Vx_E, Vy_N, Vy_S, Vy_W, Vy_E, Pa, Vxa, Vya
    # return Vx_N
end
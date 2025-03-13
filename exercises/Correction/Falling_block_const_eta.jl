using Plots
using ExtendableSparse
using GeoModBox.InitialCondition, GeoModBox.MomentumEquation.TwoD

function main()
    # Define Initial Condition ========================================== #
    # Density --- 
    #   1) block
    Ini         =   (p=:block,) 
    # ------------------------------------------------------------------- #
    # Plot Settings ===================================================== #
    Pl  =   (
        qinc    =   5,
        qsc     =   100*(60*60*24*365.25)*1e2
    )
    # ------------------------------------------------------------------- #
    # Geometry ========================================================== #
    M       =   (
        xmin    =   0.0,
        xmax    =   500.0e3,    # [ m ]
        ymin    =   -500.0e3,   # [ m ]
        ymax    =   0.0,
    )
    # -------------------------------------------------------------------- #
    # Grid =============================================================== #
    NC      =   (
        x   =   50, 
        y   =   50,
    )
    NV      =   (
        x   =   NC.x + 1,
        y   =   NC.y + 1,
    )
    Δ       =   (
        x   =   (M.xmax - M.xmin)/NC.x,
        y   =   (M.ymax - M.ymin)/NC.y,
    )
    x       =   (
        c   =   LinRange(M.xmin+Δ.x/2,M.xmax-Δ.x/2,NC.x),
        ce  =   LinRange(M.xmin - Δ.x/2.0, M.xmax + Δ.x/2.0, NC.x+2),
        v   =   LinRange(M.xmin,M.xmax,NV.x),
    )
    y       =   (
        c   =   LinRange(M.ymin+Δ.y/2,M.ymax-Δ.y/2,NC.y),
        ce  =   LinRange(M.ymin - Δ.x/2.0, M.ymax + Δ.x/2.0, NC.y+2),
        v   =   LinRange(M.ymin,M.ymax,NV.y),
    )
    x1      =   (
        c2d     =   x.c .+ 0*y.c',
        v2d     =   x.v .+ 0*y.v', 
        vx2d    =   x.v .+ 0*y.ce',
        vy2d    =   x.ce .+ 0*y.v',
    )
    x   =   merge(x,x1)
    y1      =   (
        c2d     =   0*x.c .+ y.c',
        v2d     =   0*x.v .+ y.v',
        vx2d    =   0*x.v .+ y.ce',
        vy2d    =   0*x.ce .+ y.v',
    )
    y   =   merge(y,y1)
    # -------------------------------------------------------------------- #
    # Physics ============================================================ #
    g       =   9.81

    η₀      =   1.0e21

    ρ₀      =   3200.0          #   Background density
    ρ₁      =   3300.0          #   Block density
    ρ       =   [ρ₀,ρ₁] 

    phase   =   [0,1]
    # ------------------------------------------------------------------- #
    # Allocation ======================================================== #
    D   =   (
        vx      =   zeros(Float64,NV.x,NC.y+2),
        vy      =   zeros(Float64,NC.x+2,NV.y),
        Pt      =   zeros(Float64,NC...),
        p       =   zeros(Int64,NC...),
        p_ex    =   zeros(Int64,NC.x+2,NC.y+2),
        ρ       =   zeros(Float64,NC...),
        vxc     =   zeros(Float64,NC...),
        vyc     =   zeros(Float64,NC...),
        vc      =   zeros(Float64,NC...),
    )
    # ------------------------------------------------------------------- #
    # Boundary Conditions =============================================== #
    VBC     =   (
        type    =   (E=:freeslip,W=:freeslip,S=:freeslip,N=:freeslip),
        # type    =   (E=:noslip,W=:noslip,S=:noslip,N=:noslip),
        val     =   (E=zeros(NV.y),W=zeros(NV.y),S=zeros(NV.x),N=zeros(NV.x)),
    )
    # ------------------------------------------------------------------- #
    # Initial Condition ================================================= #
    IniPhase!(Ini.p,D,M,x,y,NC;phase)
    for i in eachindex(phase)
        D.ρ[D.p.==phase[i]] .= ρ[i]
    end
    # ------------------------------------------------------------------- #
    # System of Equations =============================================== #
    # Numbering, without ghost nodes! ---
    off    = [  NV.x*NC.y,                          # vx
                NV.x*NC.y + NC.x*NV.y,              # vy
                NV.x*NC.y + NC.x*NV.y + NC.x*NC.y]  # Pt

    Num    =    (
        Vx  =   reshape(1:NV.x*NC.y, NV.x, NC.y), 
        Vy  =   reshape(off[1]+1:off[1]+NC.x*NV.y, NC.x, NV.y), 
        Pt  =   reshape(off[2]+1:off[2]+NC.x*NC.y,NC...),
    )
    rhs     =   zeros(maximum(Num.Pt))
    χ       =   zeros(maximum(Num.Pt))
    # ------------------------------------------------------------------- #
    # Assemble Coefficients ============================================= #
    K       =   Assemblyc(NC, NV, Δ, η₀, VBC, Num)
    # ------------------------------------------------------------------- #
    # Update RHS ======================================================== #
    rhs     =   updaterhsc( NC, NV, Δ, η₀, D.ρ, g, VBC, Num, rhs )
    # ------------------------------------------------------------------- #
    # Solve System of Equations ========================================= #
    χ       =   K \ rhs
    # ------------------------------------------------------------------- #
    # Update Unknown Variables ========================================== #
    D.vx[:,2:end-1]     .=  χ[Num.Vx]
    D.vy[2:end-1,:]     .=  χ[Num.Vy]
    D.Pt                .=  χ[Num.Pt]
    # ------------------------------------------------------------------- #
    # Get the velocity on the centroids ---
    for i = 1:NC.x
        for j = 1:NC.y
            D.vxc[i,j]  = (D.vx[i,j+1] + D.vx[i+1,j+1])/2
            D.vyc[i,j]  = (D.vy[i+1,j] + D.vy[i+1,j+1])/2
        end
    end
    @. D.vc        = sqrt(D.vxc^2 + D.vyc^2)

    @show(minimum(D.vc))
    @show(maximum(D.vc))

    p = heatmap(x.c./1e3,y.c./1e3,D.ρ',color=:inferno,
                    xlabel="x[km]",ylabel="y[km]",colorbar=false,
                    title="Density",
                    aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3), 
                    ylims=(M.ymin/1e3, M.ymax/1e3),
                    layout=(2,2),subplot=1)
    quiver!(p,x.c2d[1:Pl.qinc:end,1:Pl.qinc:end]./1e3,
                y.c2d[1:Pl.qinc:end,1:Pl.qinc:end]./1e3,
                quiver=(D.vx[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc,
                        D.vyc[1:Pl.qinc:end,1:Pl.qinc:end].*Pl.qsc),        
                color="white",layout=(2,2),subplot=1)
    heatmap!(p,x.c./1e3,y.c./1e3,D.vxc',
                    xlabel="x[km]",ylabel="y[km]",colorbar=false,
                    title="V_x",
                    aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),
                    ylims=(M.ymin/1e3, M.ymax/1e3),
                    layout=(2,2),subplot=3)
    heatmap!(p,x.c./1e3,y.c./1e3,D.vyc',
                    xlabel="x[km]",ylabel="y[km]",colorbar=false,
                    title="V_y",
                    aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),
                    ylims=(M.ymin/1e3, M.ymax/1e3),
                    layout=(2,2),subplot=4)
    heatmap!(p,x.c./1e3,y.c./1e3,D.Pt',
                    xlabel="x[km]",ylabel="y[km]",colorbar=false,
                    title="P_t",
                    aspect_ratio=:equal,xlims=(M.xmin/1e3, M.xmax/1e3),
                    ylims=(M.ymin/1e3, M.ymax/1e3),
                    layout=(2,2),subplot=2)
    display(p)

    savefig(p,string("./exercises/Correction/Results/09_FallingBlock_Instanteneous.png"))
end

main()
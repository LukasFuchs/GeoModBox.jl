using GeoModBox.HeatEquation
using ExactFieldSolutions
using Plots

function HeatEquation()
    # Spatial domain
    xlim = (min=-1/2, max=1/2)
    ylim = (min=-1/2, max=1/2)
    nc   = (x=100, y=100)
    nv   = (x=nc.x+1, y=nc.y+1)
    nc   = (x=nc.x+0, y=nc.y+0)
    nv   = (x=nv.x+0, y=nv.y+0)
    Δ    = (x=(xlim.max-xlim.min)/nc.x, y=(ylim.max-ylim.min)/nc.y)
    x    = (c=LinRange(xlim.min-Δ.x/2, xlim.max+Δ.x/2, nc.x), v=LinRange(xlim.min-Δ.x, xlim.max+Δ.x, nv.x))
    y    = (c=LinRange(ylim.min-Δ.y/2, ylim.max+Δ.y/2, nc.y), v=LinRange(ylim.min-Δ.y, ylim.max+Δ.y, nv.y))
    # Time domain
    nt   = 1
    Δt   = 1e-3
    t    = 0.
    # Primitive variables
    T     = zeros(nc...)
    T0    = zeros(nc...)
    Te    = zeros(nc...)
    # Derived fields
    q     = (x=zeros(nv.x, nc.x), y=zeros(nc.x, nv.x))
    # Material parameters
    ρ     = zeros(nc...)
    Cp    = zeros(nc...)
    k     = (x=zeros(nv.x, nc.x), y=zeros(nc.x, nv.x))
    # Residuals
    R     = zeros(nc...)
    # Boundary conditions
    BC   = (
        type = (W=:Neumann, E=:Neumann, S=:Neumann, N=:Neumann),
        val  = (W=zeros(nc.y), E=zeros(nc.y), S=zeros(nc.x), N=zeros(nc.x)))
    # Exact solution on cell centroids
    for i=1:nc.x, j=1:nc.y
        X = [x.c[i], y.c[j], t]
        s = Diffusion2D_Gaussian( X )
        Te[i,j] = s.u
    end
    # Exact solution on cell boundaries
    for j=1:nc.y
        X = [xlim.min, y.c[j], t]
        s = Diffusion2D_Gaussian( X )
        BC.val.W[j] = s.u
        X = [xlim.max, y.c[j], t]
        s = Diffusion2D_Gaussian( X )
        BC.val.E[j] = s.u
    end
    for i=1:nc.x
        X = [x.c[i], ylim.min, t]
        s = Diffusion2D_Gaussian( X )
        BC.val.S[i] = s.u
        X = [x.c[i], ylim.max, t]
        s = Diffusion2D_Gaussian( X )
        BC.val.N[i] = s.u
    end
    # Evaluate residual
    ComputeResiduals!(R, T, T0, ρ, Cp, k, BC, Δ, Δt)
    # Visualisation
    p1 = plot(aspect_ratio=1, xlims=(xlim...,), ylims=(ylim...,))
    p1 = heatmap!(x.c, y.c, Te')
end

HeatEquation()
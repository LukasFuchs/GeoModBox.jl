using GeoModBox.HeatEquation.TwoD
using ExactFieldSolutions, LinearAlgebra, Plots, Printf
using TimerOutputs

function HeatEquation()
    to  =   TimerOutput()
    # Spatial domain
    xlim = (min=-1/2, max=1/2)
    ylim = (min=-1/2, max=1/2)
    nc   = (x=100, y=100)
    nv   = (x=nc.x+1, y=nc.y+1)
    nc   = (x=nc.x+0, y=nc.y+0)
    nv   = (x=nv.x+0, y=nv.y+0)
    Δ    = (x=(xlim.max-xlim.min)/nc.x, y=(ylim.max-ylim.min)/nc.y)
    x    = (c=LinRange(xlim.min+Δ.x/2, xlim.max-Δ.x/2, nc.x), v=LinRange(xlim.min, xlim.max, nv.x))
    y    = (c=LinRange(ylim.min+Δ.y/2, ylim.max-Δ.y/2, nc.y), v=LinRange(ylim.min, ylim.max, nv.y))
    # Time domain
    nt   = 500
    t    = 0.
    nout = 10 
    # Iterations
    niter = 10
    ϵ     = 1e-10
    # Primitive variables
    T_ex  = zeros(nc.x+2, nc.y+2)
    T     = zeros(nc...)
    T0    = zeros(nc...)
    Te    = zeros(nc...)
    # Derived fields
    ∂T    = (∂x=zeros(nv.x, nc.x), ∂y=zeros(nc.x, nv.x))
    q     = (x=zeros(nv.x, nc.x), y=zeros(nc.x, nv.x))
    # Material parameters
    ρ     = zeros(nc...)
    Cp    = zeros(nc...)
    k     = (x=zeros(nv.x, nc.x), y=zeros(nc.x, nv.x))
    # Residuals
    R     = zeros(nc...)
    # Boundary conditions
    BC   = (
        type = (W=:Dirichlet, E=:Dirichlet, S=:Dirichlet, N=:Dirichlet),
        # type = (W=:Neumann, E=:Neumann, S=:Neumann, N=:Neumann),
        val  = (W=zeros(nc.y), E=zeros(nc.y), S=zeros(nc.x), N=zeros(nc.x)))
    # Numbering 
    Num    = (T=reshape(1:nc.x*nc.y, nc.x, nc.y),)
    # Initial conditions
    AnalyticalSolution2D!(T, x.c, y.c, t,(T0=1.0,K=1e-6,σ=0.1))
    @. k.x = 1e-6 
    @. k.y = 1e-6
    @. ρ   = 1.0
    @. Cp  = 1.0
    Δt = max(Δ...)^2/(maximum(k.x)/minimum(ρ)/minimum(Cp))/4.1
    # Time integration Loop
    @timeit to "TimeLoop" begin
    for it=1:nt
        @printf("Time step = %05d\n", it)
        t += Δt
        @. T0 = T
        # Exact solution on cell centroids
        AnalyticalSolution2D!(Te, x.c, y.c, t,(T0=1.0,K=1e-6,σ=0.1))
        # Exact solution on cell boundaries
        BoundaryConditions2D!(BC, x.c, y.c, t,(T0=1.0,K=1e-6,σ=0.1))
        # Iteration loop
        @timeit to "Iteration" begin
        for iter=1:niter
            # Evaluate residual
            @timeit to "Residual" begin
            ComputeResiduals2D!(R, T, T_ex, T0, ∂T, q, ρ, Cp, k, BC, Δ, Δt)
            @printf("||R|| = %1.4e\n", norm(R)/length(R))
            norm(R)/length(R) < ϵ ? break : nothing
            end
            # Assemble linear system
            @timeit to "Assembly" begin
            K  = AssembleMatrix2D(ρ, Cp, k, BC, Num, nc, Δ, Δt)
            end
            @timeit to "Solution" begin
            # Solve for temperature correction: Cholesky factorisation
            Kc = cholesky(K.cscmatrix)
            # Solve for temperature correction: Back substitutions
            δT = -(Kc\R[:])
            end
            # Update temperature
            @. T += δT[Num.T]
        end
        end
        # Visualisation
        if mod(it, nout)==0
            p1 = plot(aspect_ratio=1, xlims=(xlim...,), ylims=(ylim...,))
            p1 = heatmap!(x.c, y.c, Te', title="Analytics")
            p2 = plot(aspect_ratio=1, xlims=(xlim...,), ylims=(ylim...,))
            p2 = heatmap!(x.c, y.c, T', title="Numerics")
            p3 = plot(aspect_ratio=1, xlims=(xlim...,), ylims=(ylim...,))
            p3 = heatmap!(x.c, y.c, (abs.(T-Te))', title="Error")
            display(plot(p1, p2, p3, layout=(2,2)))
        end
    end
    end
    display(to)
end

HeatEquation()
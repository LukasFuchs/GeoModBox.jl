using GeoModBox.HeatEquation.TwoD
using Plots
using TimerOutputs
using ExactFieldSolutions

function HeatEquation()

    to = TimerOutput()

    # -------------------------------------------------------------------------
    # Spatial domain
    # -------------------------------------------------------------------------

    xlim = (min = -1 / 2, max = 1 / 2)
    ylim = (min = -1 / 2, max = 1 / 2)

    nc = (x = 100, y = 100)
    nv = (x = nc.x + 1, y = nc.y + 1)

    Δ = (
        x = (xlim.max - xlim.min) / nc.x,
        y = (ylim.max - ylim.min) / nc.y
    )

    x = (
        c = LinRange(xlim.min + Δ.x / 2, xlim.max - Δ.x / 2, nc.x),
        v = LinRange(xlim.min, xlim.max, nv.x)
    )

    y = (
        c = LinRange(ylim.min + Δ.y / 2, ylim.max - Δ.y / 2, nc.y),
        v = LinRange(ylim.min, ylim.max, nv.y)
    )

    # -------------------------------------------------------------------------
    # Time domain
    # -------------------------------------------------------------------------

    nt   = 400
    t    = 0.0
    nout = 10

    # -------------------------------------------------------------------------
    # Primitive variables
    # -------------------------------------------------------------------------

    T_ex = zeros(nc.x + 2, nc.y + 2)
    T    = zeros(nc...)
    Te   = zeros(nc...)

    # -------------------------------------------------------------------------
    # Derived fields
    # -------------------------------------------------------------------------

    ∂T = (
        ∂x = zeros(nv.x, nc.y),
        ∂y = zeros(nc.x, nv.y)
    )

    q = (
        x = zeros(nv.x, nc.y),
        y = zeros(nc.x, nv.y)
    )

    # -------------------------------------------------------------------------
    # Material parameters
    # -------------------------------------------------------------------------

    ρ  = zeros(nc...)
    Cp = zeros(nc...)

    k = (
        x = zeros(nv.x, nc.y),
        y = zeros(nc.x, nv.y)
    )

    # -------------------------------------------------------------------------
    # Boundary conditions
    # -------------------------------------------------------------------------

    BC = (
        type = (
            W = :Dirichlet,
            E = :Dirichlet,
            S = :Dirichlet,
            N = :Dirichlet
        ),
        val = (
            W = zeros(nc.y),
            E = zeros(nc.y),
            S = zeros(nc.x),
            N = zeros(nc.x)
        )
    )

    # -------------------------------------------------------------------------
    # Initial conditions
    # -------------------------------------------------------------------------

    AnalyticalSolution2D!(T, x.c, y.c, t, (T0 = 1.0, K = 1e-6, σ = 0.1))

    @. k.x = 1e-6
    @. k.y = 1e-6
    @. ρ   = 1.0
    @. Cp  = 1.0

    κmax = max(maximum(k.x), maximum(k.y)) /
           (minimum(ρ) * minimum(Cp))

    Δt = 1 / (2.05 * κmax * (1 / Δ.x^2 + 1 / Δ.y^2))

    # -------------------------------------------------------------------------
    # Time integration
    # -------------------------------------------------------------------------

    @timeit to "TimeLoop" begin

        @views for it = 1:nt

            # -----------------------------------------------------------------
            # Exact solution on the model boundaries
            # -----------------------------------------------------------------

            BoundaryConditions2D!(BC, x.c, y.c, t, (T0 = 1.0, K = 1e-6, σ = 0.1))

            # -----------------------------------------------------------------
            # Copy centroid temperatures to expanded grid
            # -----------------------------------------------------------------

            @. T_ex[2:end-1, 2:end-1] = T

            # -----------------------------------------------------------------
            # Update ghost nodes
            # -----------------------------------------------------------------

            @. T_ex[1, 2:end-1] =
                (BC.type.W == :Dirichlet) *
                (2 * BC.val.W - T_ex[2, 2:end-1])

            @. T_ex[end, 2:end-1] =
                (BC.type.E == :Dirichlet) *
                (2 * BC.val.E - T_ex[end-1, 2:end-1])

            @. T_ex[2:end-1, 1] =
                (BC.type.S == :Dirichlet) *
                (2 * BC.val.S - T_ex[2:end-1, 2])

            @. T_ex[2:end-1, end] =
                (BC.type.N == :Dirichlet) *
                (2 * BC.val.N - T_ex[2:end-1, end-1])

            # -----------------------------------------------------------------
            # Temperature gradients
            # -----------------------------------------------------------------

            @. ∂T.∂x =
                (T_ex[2:end, 2:end-1] -
                 T_ex[1:end-1, 2:end-1]) / Δ.x

            @. ∂T.∂y =
                (T_ex[2:end-1, 2:end] -
                 T_ex[2:end-1, 1:end-1]) / Δ.y

            # -----------------------------------------------------------------
            # Conductive heat fluxes
            # -----------------------------------------------------------------

            @. q.x = -k.x * ∂T.∂x
            @. q.y = -k.y * ∂T.∂y

            # -----------------------------------------------------------------
            # Forward Euler update
            # -----------------------------------------------------------------

            @. T -= Δt / (ρ * Cp) * (
                (q.x[2:end, :] - q.x[1:end-1, :]) / Δ.x +
                (q.y[:, 2:end] - q.y[:, 1:end-1]) / Δ.y
            )

            # Advance time

            t += Δt

            # -----------------------------------------------------------------
            # Exact solution on cell centroids
            # -----------------------------------------------------------------

            AnalyticalSolution2D!(Te, x.c, y.c, t, (T0 = 1.0, K = 1e-6, σ = 0.1))

            # -----------------------------------------------------------------
            # Visualisation
            # -----------------------------------------------------------------

            if mod(it, nout) == 0

                p1 = plot(aspect_ratio = 1, xlims = (xlim...,), ylims = (ylim...,))
                p1 = heatmap!(x.c, y.c, Te', title = "Analytical")

                p2 = plot(aspect_ratio = 1, xlims = (xlim...,), ylims = (ylim...,))
                p2 = heatmap!(x.c, y.c, T', title = "Numerics")

                p3 = plot(aspect_ratio = 1, xlims = (xlim...,), ylims = (ylim...,))
                p3 = heatmap!(x.c, y.c, abs.(T .- Te)', title = "Error")

                display(plot(p1, p2, p3, layout = (2, 2)))

            end
        end
    end

    display(to)

end

HeatEquation()
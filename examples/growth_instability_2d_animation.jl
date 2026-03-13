using Printf
using CairoMakie
using CartesianGrids
using PenguinBCs
using PenguinExtendedStefan

# 2D growth-instability animation: temperature + concentration.
# Output: MP4 animation rendered with CairoMakie.

function build_problem(; nx=65, ny=65)
    Lx = 2.0
    Ly = 1.0
    xΓ0 = 0.0
    amp = 0.03

    grid = CartesianGrid((-Lx / 2, -Ly / 2), (Lx / 2, Ly / 2), (nx, ny))
    rep = LevelSetRep(grid, (x, y) -> x - (xΓ0 + amp * sin(2 * pi * y / Ly)))

    bcT = BorderConditions(; left=Dirichlet(1.0), right=Dirichlet(-0.3), bottom=Neumann(0.0), top=Neumann(0.0))
    bcC = BorderConditions(; left=Dirichlet(1.0), right=Dirichlet(0.0), bottom=Neumann(0.0), top=Neumann(0.0))

    params = ExtendedStefanParams(
        1.0,
        1.0,
        1.5,
        AlloyEquilibrium(0.7, 0.0, -0.25);
        kappa_l=1.0,
        kappa_s=0.9,
        D_l=0.1,
        D_s=0.03,
    )
    opts = ExtendedStefanOptions(
        ;
        scheme=:BE,
        reinit=true,
        reinit_every=5,
        coupling_max_iter=12,
        coupling_tol=1e-7,
        coupling_reltol=1e-7,
        coupling_damping=0.6,
    )

    prob = ExtendedStefanProblem(grid, bcT, bcC, params, rep, opts)

    solver = build_solver(
        prob;
        Tlω0=(x, y, t) -> 0.5 - 0.8 * x,
        Tlγ0=0.0,
        Tsω0=(x, y, t) -> -0.2 - 0.2 * x,
        Tsγ0=0.0,
        Clω0=(x, y, t) -> 0.5 - 0.4 * x,
        Clγ0=0.2,
        Csω0=(x, y, t) -> 0.1 - 0.1 * x,
        Csγ0=0.14,
        speed0=0.0,
    )

    return solver
end

function run_with_retry!(solver; t0=0.0, tf=0.03, dt0=0.005)
    dt = dt0
    success = false
    last_err = nothing
    out = nothing
    while dt >= 5e-4 && !success
        try
            out = solve!(solver, (t0, tf); dt=dt, save_history=true)
            success = all(isfinite, out.state.Tlω) && all(isfinite, out.state.Clω)
        catch err
            last_err = err
            success = false
        end
        success || (dt *= 0.5)
    end
    success || error("Growth-instability run failed (last error: $(last_err))")
    return out, dt
end

function interface_amplitude(phi::AbstractMatrix{T}, x::AbstractVector{T}) where {T}
    n1, n2 = size(phi)
    xs = T[]
    @inbounds for j in 1:n2
        for i in 1:(n1 - 1)
            a = phi[i, j]
            b = phi[i + 1, j]
            if isfinite(a) && isfinite(b) && a * b <= 0
                θ = (abs(a) + abs(b)) > 0 ? abs(a) / (abs(a) + abs(b)) : convert(T, 0.5)
                push!(xs, (one(T) - θ) * x[i] + θ * x[i + 1])
                break
            end
        end
    end
    isempty(xs) && return T(NaN)
    return convert(T, 0.5) * (maximum(xs) - minimum(xs))
end

function composite_fields(frame, nshape)
    n1, n2 = nshape
    ϕ = reshape(frame.phi, n1, n2)

    Tl = reshape(frame.Tlω, n1, n2)
    Ts = reshape(frame.Tsω, n1, n2)
    Cl = reshape(frame.Clω, n1, n2)
    Cs = reshape(frame.Csω, n1, n2)

    T = similar(Tl)
    C = similar(Cl)
    @inbounds for j in 1:n2, i in 1:n1
        if ϕ[i, j] <= 0
            T[i, j] = Tl[i, j]
            C[i, j] = Cl[i, j]
        else
            T[i, j] = Ts[i, j]
            C[i, j] = Cs[i, j]
        end
    end
    return T, C, ϕ
end

function collect_movie_fields(history, nshape)
    Tframes = Matrix{Float64}[]
    Cframes = Matrix{Float64}[]
    Pframes = Matrix{Float64}[]
    times = Float64[]

    for fr in history
        T, C, P = composite_fields(fr, nshape)
        push!(Tframes, Matrix{Float64}(T))
        push!(Cframes, Matrix{Float64}(C))
        push!(Pframes, Matrix{Float64}(P))
        push!(times, fr.t)
    end

    return Tframes, Cframes, Pframes, times
end

function finite_extrema(frames)
    fmin = Inf
    fmax = -Inf
    for A in frames
        @inbounds for v in A
            isfinite(v) || continue
            fmin = min(fmin, v)
            fmax = max(fmax, v)
        end
    end
    isfinite(fmin) && isfinite(fmax) || error("No finite values found for color range")
    return fmin, fmax
end

function animate_growth_instability(; output="growth_instability_2d_temperature_concentration.mp4", framerate=20)
    solver = build_problem()
    out, dt = run_with_retry!(solver)

    grid = solver.problem.grid
    nshape = grid.n
    x = collect(CartesianGrids.grid1d(grid, 1))
    y = collect(CartesianGrids.grid1d(grid, 2))

    Tframes, Cframes, Pframes, times = collect_movie_fields(out.history, nshape)
    nframes = length(times)

    amps = [interface_amplitude(Pframes[i], x) for i in 1:nframes]
    a0 = amps[1]
    g = amps[end] / a0

    tmin, tmax = finite_extrema(Tframes)
    cmin, cmax = finite_extrema(Cframes)

    Tobs = Observable(Tframes[1])
    Cobs = Observable(Cframes[1])
    Pobs = Observable(Pframes[1])
    title_obs = Observable(@sprintf("t = %.4f", times[1]))

    fig = Figure(size=(1300, 620), fontsize=18)
    axT = Axis(fig[1, 1], title="Temperature", xlabel="x", ylabel="y")
    axC = Axis(fig[1, 2], title="Concentration", xlabel="x", ylabel="y")

    hmT = heatmap!(axT, x, y, Tobs; colormap=:thermal, colorrange=(tmin, tmax), nan_color=:transparent)
    hmC = heatmap!(axC, x, y, Cobs; colormap=:haline, colorrange=(cmin, cmax), nan_color=:transparent)

    contour!(axT, x, y, Pobs; levels=[0.0], color=:black, linewidth=2)
    contour!(axC, x, y, Pobs; levels=[0.0], color=:black, linewidth=2)

    Colorbar(fig[2, 1], hmT, vertical=false, label="T")
    Colorbar(fig[2, 2], hmC, vertical=false, label="C")
    Label(fig[0, 1:2], title_obs, fontsize=22)

    record(fig, output, 1:nframes; framerate=framerate) do i
        Tobs[] = Tframes[i]
        Cobs[] = Cframes[i]
        Pobs[] = Pframes[i]
        title_obs[] = @sprintf("Growth instability: t = %.4f", times[i])
    end

    metrics = extended_stefan_residual_metrics(out.solver)
    println("Saved animation: $(abspath(output))")
    println("  frames   = ", nframes)
    println("  dt used  = ", dt)
    println("  amp      = ", a0, " -> ", amps[end], " (growth = ", g, ")")
    println("  max speed= ", metrics.max_speed)
end

output = isempty(ARGS) ? "growth_instability_2d_temperature_concentration.mp4" : ARGS[1]
animate_growth_instability(; output=output)

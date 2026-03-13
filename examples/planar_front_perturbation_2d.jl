using CartesianGrids
using PenguinBCs
using PenguinExtendedStefan

# Exploratory example: planar front with a small sinusoidal perturbation.
# This is not a validated dendrite benchmark.

Lx = 2.0
Ly = 1.0
xΓ0 = 0.0
amp = 0.03

grid = CartesianGrid((-Lx / 2, -Ly / 2), (Lx / 2, Ly / 2), (65, 65))
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

function interface_amplitude(phi, grid)
    n1, n2 = grid.n
    x = collect(CartesianGrids.grid1d(grid, 1))
    xs = Float64[]
    @inbounds for j in 1:n2
        for i in 1:(n1 - 1)
            a = phi[i, j]
            b = phi[i + 1, j]
            if isfinite(a) && isfinite(b) && a * b <= 0
                θ = (abs(a) + abs(b)) > 0 ? abs(a) / (abs(a) + abs(b)) : 0.5
                push!(xs, (1 - θ) * x[i] + θ * x[i + 1])
                break
            end
        end
    end
    isempty(xs) && return NaN
    return 0.5 * (maximum(xs) - minimum(xs))
end

function run_with_retry!(solver)
    dt = 0.005
    success = false
    last_err = nothing
    out = nothing
    while dt >= 5e-4 && !success
        try
            out = solve!(solver, (0.0, 0.03); dt=dt, save_history=true)
            success = all(isfinite, out.state.Tlω) && all(isfinite, out.state.Clω)
        catch err
            last_err = err
            success = false
        end
        success || (dt *= 0.5)
    end
    success || error("Perturbation example failed to converge (last error: $(last_err))")
    return out, dt
end

out, dt = run_with_retry!(solver)
amp0 = interface_amplitude(reshape(out.history[1].phi, grid.n...), grid)
amp1 = interface_amplitude(reshape(out.history[end].phi, grid.n...), grid)

metrics = extended_stefan_residual_metrics(out.solver)

println("Planar-front perturbation (2D exploratory)")
println("  final time              = ", out.state.t)
println("  used dt                 = ", dt)
println("  interface amplitude     = ", amp0, " -> ", amp1, " (growth = ", amp1 / amp0, ")")
println("  max speed               = ", metrics.max_speed)
println("  nonlinear iterations    = ", metrics.nonlinear_iteration_count)
println("  thermal residual (last) = ", metrics.thermal_stefan_residual_norm)
println("  solute residual (last)  = ", metrics.solute_residual_norm)

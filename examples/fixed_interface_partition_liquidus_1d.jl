using CartesianGrids
using PenguinBCs
using PenguinExtendedStefan

grid = CartesianGrid((0.0,), (1.0,), (129,))
rep = LevelSetRep(grid, x -> x - 0.53)

bcT = BorderConditions(; left=Dirichlet(0.8), right=Dirichlet(0.2))
bcC = BorderConditions(; left=Dirichlet(0.6), right=Dirichlet(0.0))

k = 0.65
Tm = 0.0
mliq = -0.25

params = ExtendedStefanParams(
    1.0,
    1.0,
    1e12,
    AlloyEquilibrium(k, Tm, mliq);
    kappa_l=1.0,
    kappa_s=1.0,
    D_l=0.1,
    D_s=0.05,
)
opts = ExtendedStefanOptions(; scheme=:BE, reinit=false, coupling_max_iter=12, coupling_damping=0.7)
prob = ExtendedStefanProblem(grid, bcT, bcC, params, rep, opts)

solver = build_solver(
    prob;
    Tlω0=0.5,
    Tlγ0=0.4,
    Tsω0=0.4,
    Tsγ0=0.4,
    Clω0=0.3,
    Clγ0=0.25,
    Csω0=0.15,
    Csγ0=0.18,
    speed0=0.0,
)

solve!(solver, (0.0, 0.02); dt=0.01, save_history=false)

alg = interface_algebra_residuals(solver)
println("Fixed-interface partition/liquidus check (1D)")
println("  max |TlΓ - TsΓ|             = ", alg.continuity.maxabs)
println("  max |CsΓ - k ClΓ|           = ", alg.partition.maxabs)
println("  max |TlΓ - (Tm + m ClΓ)|    = ", alg.liquidus.maxabs)
println("  max |speed|                 = ", maximum(abs.(vec(solver.state.speed_full))))

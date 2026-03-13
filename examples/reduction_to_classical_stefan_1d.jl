using CartesianGrids
using PenguinBCs
using PenguinExtendedStefan

function run_case(c0)
    grid = CartesianGrid((0.0,), (1.0,), (129,))
    rep = LevelSetRep(grid, x -> x - 0.46)

    bcT = BorderConditions(; left=Dirichlet(1.0), right=Dirichlet(-0.2))
    bcC = BorderConditions(; left=Dirichlet(c0), right=Dirichlet(c0))

    params = ExtendedStefanParams(
        1.0,
        1.0,
        8.0,
        AlloyEquilibrium(1.0, 0.0, 0.0);
        kappa_l=1.0,
        kappa_s=1.0,
        D_l=0.1,
        D_s=0.1,
        source_T_l=0.0,
        source_T_s=0.0,
        source_C_l=0.0,
        source_C_s=0.0,
    )
    opts = ExtendedStefanOptions(; scheme=:BE, reinit=false, coupling_max_iter=10)
    prob = ExtendedStefanProblem(grid, bcT, bcC, params, rep, opts)

    solver = build_solver(
        prob;
        Tlω0=(x, t) -> 1 - x,
        Tlγ0=0.0,
        Tsω0=(x, t) -> -0.2 * (x - 0.46),
        Tsγ0=0.0,
        Clω0=c0,
        Clγ0=c0,
        Csω0=c0,
        Csγ0=c0,
        speed0=0.0,
    )

    solve!(solver, (0.0, 0.05); dt=0.01, save_history=false)
    return solver
end

s0 = run_case(0.0)
s1 = run_case(0.8)

println("Reduction-to-classical-Stefan check (1D)")
println("  reduction_to_thermal_stefan = ", PenguinExtendedStefan.reduction_to_thermal_stefan(s0.problem.params))
println("  max |Cl(c0=0)-0|            = ", maximum(abs.(s0.state.Clω)))
println("  max |Cl(c0=0.8)-0.8|        = ", maximum(abs.(s1.state.Clω .- 0.8)))
println("  max |Tl(c0=0)-Tl(c0=0.8)|   = ", maximum(abs.(s0.state.Tlω .- s1.state.Tlω)))
println("  max |speed diff|            = ", maximum(abs.(vec(s0.state.speed_full .- s1.state.speed_full))))

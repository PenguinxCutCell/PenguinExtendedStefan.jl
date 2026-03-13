function _run_decoupled_case(c0)
    grid = CartesianGrid((0.0,), (1.0,), (81,))
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
    opts = ExtendedStefanOptions(; scheme=:BE, reinit=false, coupling_max_iter=10, coupling_tol=1e-8, coupling_reltol=1e-8)
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

    solve!(solver, (0.0, 0.04); dt=0.01, save_history=false)
    return solver
end

@testset "Reduction to classical Stefan behavior (1D)" begin
    solver0 = _run_decoupled_case(0.0)
    solver1 = _run_decoupled_case(0.8)

    p = solver0.problem.params
    @test PenguinExtendedStefan.reduction_to_thermal_stefan(p)

    @test all(isfinite, solver0.state.Clω)
    @test all(isfinite, solver0.state.Csω)
    @test all(isfinite, solver1.state.Clω)
    @test all(isfinite, solver1.state.Csω)

    @test maximum(abs.(solver0.state.Tlω .- solver1.state.Tlω)) < 1e-9
    @test maximum(abs.(solver0.state.Tsω .- solver1.state.Tsω)) < 1e-9
    @test maximum(abs.(vec(solver0.state.speed_full .- solver1.state.speed_full))) < 1e-9

    alg0 = interface_algebra_residuals(solver0)
    @test alg0.liquidus.maxabs < 1e-8
    @test alg0.partition.maxabs < 1e-8

    metrics = extended_stefan_residual_metrics(solver0)
    @test metrics.max_speed > 0.0
end

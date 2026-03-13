@testset "Curved interface physics check (2D smoke)" begin
    grid = CartesianGrid((-1.0, -1.0), (1.0, 1.0), (33, 33))
    rep = LevelSetRep(grid, (x, y) -> hypot(x, y) - 0.45)

    bcT = BorderConditions(; left=Dirichlet(1.0), right=Dirichlet(-0.2), bottom=Dirichlet(0.8), top=Dirichlet(-0.1))
    bcC = BorderConditions(; left=Dirichlet(0.7), right=Dirichlet(0.1), bottom=Dirichlet(0.6), top=Dirichlet(0.05))

    params = ExtendedStefanParams(
        1.0,
        1.0,
        1e12,
        AlloyEquilibrium(0.65, 0.0, -0.25);
        kappa_l=1.0,
        kappa_s=0.9,
        D_l=0.1,
        D_s=0.05,
    )
    opts = ExtendedStefanOptions(; scheme=:BE, reinit=false, coupling_max_iter=8, coupling_tol=1e-8, coupling_reltol=1e-8)
    prob = ExtendedStefanProblem(grid, bcT, bcC, params, rep, opts)

    solver = build_solver(
        prob;
        Tlω0=0.2,
        Tlγ0=0.1,
        Tsω0=-0.1,
        Tsγ0=0.1,
        Clω0=0.3,
        Clγ0=0.2,
        Csω0=0.1,
        Csγ0=0.13,
        speed0=0.0,
    )

    step!(solver, 0.005)

    @test all(isfinite, solver.state.Tlω)
    @test all(isfinite, solver.state.Tsω)
    @test all(isfinite, solver.state.Clω)
    @test all(isfinite, solver.state.Csω)

    @test maximum(abs.(solver.state.Tlω)) < 10.0
    @test maximum(abs.(solver.state.Clω)) < 10.0

    alg = interface_algebra_residuals(solver)
    @test alg.continuity.maxabs < 1e-7
    @test alg.partition.maxabs < 1e-7
    @test alg.liquidus.maxabs < 1e-7

    m = extended_stefan_residual_metrics(solver)
    @test isfinite(m.solute_residual_norm)
    @test isfinite(m.thermal_stefan_residual_norm)
end

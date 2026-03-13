@testset "Prescribed circular motion (2D smoke)" begin
    case = prescribed_circle_motion_case(; center=(0.02, -0.01), R0=0.32, V0=0.03)
    t0 = 0.0
    tf = 0.04
    dt = 0.004

    grid = CartesianGrid((-1.0, -1.0), (1.0, 1.0), (33, 33))
    rep = LevelSetRep(grid, case.phi0)

    bcT = BorderConditions(
        ;
        left=Dirichlet((x, y, t) -> case.Tl_exact(x, y, t)),
        right=Dirichlet((x, y, t) -> case.Ts_exact(x, y, t)),
        bottom=Dirichlet((x, y, t) -> case.Tl_exact(x, y, t)),
        top=Dirichlet((x, y, t) -> case.Ts_exact(x, y, t)),
    )
    bcC = BorderConditions(
        ;
        left=Dirichlet((x, y, t) -> case.Cl_exact(x, y, t)),
        right=Dirichlet((x, y, t) -> case.Cs_exact(x, y, t)),
        bottom=Dirichlet((x, y, t) -> case.Cl_exact(x, y, t)),
        top=Dirichlet((x, y, t) -> case.Cs_exact(x, y, t)),
    )

    params = ExtendedStefanParams(
        1.0,
        1.0,
        1.0,
        AlloyEquilibrium(0.75, 0.0, -0.2);
        kappa_l=1.0,
        kappa_s=0.9,
        D_l=0.1,
        D_s=0.04,
        source_T_l=case.source_T_l,
        source_T_s=case.source_T_s,
        source_C_l=case.source_C_l,
        source_C_s=case.source_C_s,
    )
    opts = ExtendedStefanOptions(; scheme=:CN, reinit=false, coupling_max_iter=2, coupling_damping=1.0)
    prob = ExtendedStefanProblem(grid, bcT, bcC, params, rep, opts)

    solver = build_solver(
        prob;
        t0=t0,
        Tlω0=(x, y, t) -> case.Tl_exact(x, y, t0),
        Tlγ0=(x, y, t) -> case.Tl_exact(x, y, t0),
        Tsω0=(x, y, t) -> case.Ts_exact(x, y, t0),
        Tsγ0=(x, y, t) -> case.Ts_exact(x, y, t0),
        Clω0=(x, y, t) -> case.Cl_exact(x, y, t0),
        Clγ0=(x, y, t) -> case.Cl_exact(x, y, t0),
        Csω0=(x, y, t) -> case.Cs_exact(x, y, t0),
        Csγ0=(x, y, t) -> case.Cs_exact(x, y, t0),
        speed0=0.0,
    )

    solve_prescribed_motion!(solver, (t0, tf); dt=dt, speed_callback=case.speed, save_history=false)

    @test all(isfinite, solver.state.Tlω)
    @test all(isfinite, solver.state.Tsω)
    @test all(isfinite, solver.state.Clω)
    @test all(isfinite, solver.state.Csω)

    vmax = maximum(abs.(vec(solver.state.speed_full)))
    @test isapprox(vmax, case.V(tf); atol=1e-12)

    alg = interface_algebra_residuals(solver)
    @test alg.continuity.maxabs < 1e-7
    @test alg.partition.maxabs < 1e-7
    @test alg.liquidus.maxabs < 1e-7
end

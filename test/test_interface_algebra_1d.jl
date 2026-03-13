@testset "Fixed-interface alloy algebra (1D)" begin
    grid = CartesianGrid((0.0,), (1.0,), (97,))
    rep = LevelSetRep(grid, x -> x - 0.53)

    bcT = BorderConditions(; left=Dirichlet(0.85), right=Dirichlet(0.15))
    bcC = BorderConditions(; left=Dirichlet(0.70), right=Dirichlet(0.05))

    k = 0.62
    Tm = 0.0
    mliq = -0.30
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
    opts = ExtendedStefanOptions(; scheme=:BE, reinit=false, coupling_max_iter=12, coupling_tol=1e-8, coupling_reltol=1e-8)
    prob = ExtendedStefanProblem(grid, bcT, bcC, params, rep, opts)

    solver = build_solver(
        prob;
        Tlω0=0.5,
        Tlγ0=0.4,
        Tsω0=0.3,
        Tsγ0=0.4,
        Clω0=0.30,
        Clγ0=0.25,
        Csω0=0.18,
        Csγ0=0.15,
        speed0=0.0,
    )

    solve!(solver, (0.0, 0.02); dt=0.01, save_history=false)

    alg = interface_algebra_residuals(solver)
    @test alg.continuity.count > 0
    @test alg.continuity.maxabs < 1e-8
    @test alg.partition.maxabs < 1e-8
    @test alg.liquidus.maxabs < 1e-8

    st = solver.state
    @test maximum(abs.(st.Tlγ .- st.Tsγ)) < 1e-8
    @test maximum(abs.(st.Csγ .- k .* st.Clγ)) < 1e-8
    @test maximum(abs.(st.Tlγ .- (Tm .+ mliq .* st.Clγ))) < 1e-8
    @test maximum(abs.(vec(st.speed_full))) < 1e-4

    m = extended_stefan_residual_metrics(solver)
    @test isfinite(m.thermal_stefan_residual_norm)
    @test isfinite(m.solute_residual_norm)
end

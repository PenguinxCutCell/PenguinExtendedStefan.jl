function _run_mk_level(n)
    t0 = 0.05
    tf = 0.08
    L = 0.8

    solver, pars = mk_reference_problem(n; t0=t0, L=L)
    dx = 2L / (n - 1)
    dt = 0.2 * dx

    solve!(solver, (t0, tf); dt=dt, save_history=false)

    prob = solver.problem
    h_num = PenguinExtendedStefan.interface_position_from_levelset_1d(prob.interface_rep, prob.grid)
    h_ref = martin_kauffman_interface_position(tf, pars)

    capTl = something(solver.cache.modelT.cap1_slab)
    capCl = something(solver.cache.modelC.cap1_slab)

    eT_l = profile_error_on_phase(solver.state.Tlω, (x, t) -> martin_kauffman_fields(x, tf, pars).T, capTl, tf)
    eC_l = profile_error_on_phase(solver.state.Clω, (x, t) -> martin_kauffman_fields(x, tf, pars).C, capCl, tf)

    return (
        dx=dx,
        e_h=abs(h_num - h_ref),
        eT=eT_l.L2,
        eC=eC_l.L2,
        solver=solver,
    )
end

@testset "Martin-Kauffman 1977 benchmark (1D)" begin
    pars = martin_kauffman_reference()
    @test isapprox(pars.S, 2.5; atol=1e-12)
    @test isapprox(pars.tau, 0.1; atol=1e-12)
    @test isapprox(pars.Lambda, 0.4; atol=1e-12)
    @test isapprox(pars.A, 0.90954; atol=1e-8)
    @test isapprox(pars.B, 0.44748; atol=1e-8)
    @test isapprox(pars.alpha, 0.19742; atol=1e-8)

    tchk = 0.08
    @test isapprox(
        martin_kauffman_interface_position(tchk, pars),
        2 * pars.alpha * sqrt(pars.kappaT * tchk);
        atol=1e-14,
    )

    f = martin_kauffman_fields(0.0, tchk, pars)
    @test isfinite(f.T)
    @test isfinite(f.C)
    @test isfinite(f.h)

    out49 = _run_mk_level(49)
    out65 = _run_mk_level(65)
    out97 = _run_mk_level(97)

    @test isfinite(out49.e_h) && isfinite(out65.e_h) && isfinite(out97.e_h)
    @test out65.e_h < out49.e_h
    @test out97.e_h < out65.e_h

    p12 = log(out49.e_h / out65.e_h) / log(out49.dx / out65.dx)
    p23 = log(out65.e_h / out97.e_h) / log(out65.dx / out97.dx)
    @test p12 > 0.9
    @test p23 > 0.9

    @test out97.e_h < 5e-4
    @test out97.eT < 1e-2
    @test out97.eC < 1e-2

    metrics = extended_stefan_residual_metrics(out97.solver)
    @test metrics.nonlinear_iteration_count <= out97.solver.problem.options.coupling_max_iter
end

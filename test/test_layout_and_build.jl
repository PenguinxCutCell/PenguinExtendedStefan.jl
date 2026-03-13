@testset "Layout and build API" begin
    alloy = AlloyEquilibrium(0.7, 0.0, -0.2)
    params = ExtendedStefanParams(1.0, 1.0, 5.0, alloy; kappa_l=1.2, kappa_s=0.9, D_l=0.11, D_s=0.05)
    opts = ExtendedStefanOptions(; scheme=:BE, coupling_max_iter=6)
    @test params.rho_cp_l == 1.0
    @test params.rho_cp_s == 1.0
    @test opts.scheme === :BE
    @test opts.coupling_max_iter == 6

    prob = basic_problem_1d(n=41)
    @test prob isa ExtendedStefanProblem

    solver = make_solver_from_problem(prob)
    nt = prod(prob.grid.n)

    @test length(solver.state.Tlω) == nt
    @test length(solver.state.Tlγ) == nt
    @test length(solver.state.Tsω) == nt
    @test length(solver.state.Tsγ) == nt
    @test length(solver.state.Clω) == nt
    @test length(solver.state.Clγ) == nt
    @test length(solver.state.Csω) == nt
    @test length(solver.state.Csγ) == nt
    @test size(solver.state.speed_full) == prob.grid.n
    @test length(solver.state.logs[:times]) == 1

    lay = solver.cache.layout
    idx = Int[]
    append!(idx, lay.Tlω)
    append!(idx, lay.Tlγ)
    append!(idx, lay.Tsω)
    append!(idx, lay.Tsγ)
    append!(idx, lay.Clω)
    append!(idx, lay.Clγ)
    append!(idx, lay.Csω)
    append!(idx, lay.Csγ)
    @test length(idx) == 8 * nt
    @test length(unique(idx)) == 8 * nt
    @test minimum(idx) == 1
    @test maximum(idx) == 8 * nt

    t = solver.state.t
    dt = 0.01
    phi = PenguinExtendedStefan.phi_values(prob.interface_rep)
    sys = PenguinExtendedStefan.assemble_monolithic_system!(
        solver.cache,
        phi,
        phi,
        solver.state.speed_full,
        solver.state.Tlω,
        solver.state.Tlγ,
        solver.state.Tsω,
        solver.state.Tsγ,
        solver.state.Clω,
        solver.state.Clγ,
        solver.state.Csω,
        solver.state.Csγ,
        t,
        dt,
        prob,
    )
    @test size(sys.A) == (8 * nt, 8 * nt)
    @test length(sys.b) == 8 * nt
    @test all(isfinite, sys.A.nzval)
    @test all(isfinite, sys.b)

    t0 = solver.state.t
    step!(solver, dt)
    @test solver.state.t == t0 + dt
    @test all(isfinite, solver.state.Tlω)
    @test all(isfinite, solver.state.Clω)
    @test !isempty(solver.state.logs[:speed_max])
    @test !isempty(solver.state.logs[:nonlinear_iters])

    out = solve!(solver, (solver.state.t, solver.state.t + 0.02); dt=0.01, save_history=true)
    @test isfinite(out.state.t)
    @test length(out.history) == 3
    @test all(isfinite, out.state.Tsω)
    @test all(isfinite, out.state.Csω)
end

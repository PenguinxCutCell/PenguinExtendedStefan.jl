function _run_mms_fixed_1d(n)
    case = manufactured_fixed_interface_1d_case(; xΓ=0.41)
    t0 = 0.2
    tf = 0.32
    h = 1 / (n - 1)
    dt = 0.4 * h

    grid = CartesianGrid((0.0,), (1.0,), (n,))
    rep = LevelSetRep(grid, case.phi0)

    bcT = BorderConditions(; left=Dirichlet((x, t) -> case.Tl_exact(x, t)), right=Dirichlet((x, t) -> case.Ts_exact(x, t)))
    bcC = BorderConditions(; left=Dirichlet((x, t) -> case.Cl_exact(x, t)), right=Dirichlet((x, t) -> case.Cs_exact(x, t)))

    params = ExtendedStefanParams(
        1.0,
        1.0,
        1e12,
        AlloyEquilibrium(case.k, case.Tm, case.m_liquidus);
        kappa_l=1.0,
        kappa_s=0.8,
        D_l=0.12,
        D_s=0.05,
        source_T_l=case.source_T_l,
        source_T_s=case.source_T_s,
        source_C_l=case.source_C_l,
        source_C_s=case.source_C_s,
    )
    opts = ExtendedStefanOptions(; scheme=:CN, reinit=false, coupling_max_iter=10, coupling_tol=1e-8, coupling_reltol=1e-8)
    prob = ExtendedStefanProblem(grid, bcT, bcC, params, rep, opts)

    solver = build_solver(
        prob;
        t0=t0,
        Tlω0=(x, t) -> case.Tl_exact(x, t0),
        Tlγ0=(x, t) -> case.Tl_exact(x, t0),
        Tsω0=(x, t) -> case.Ts_exact(x, t0),
        Tsγ0=(x, t) -> case.Ts_exact(x, t0),
        Clω0=(x, t) -> case.Cl_exact(x, t0),
        Clγ0=(x, t) -> case.Cl_exact(x, t0),
        Csω0=(x, t) -> case.Cs_exact(x, t0),
        Csγ0=(x, t) -> case.Cs_exact(x, t0),
        speed0=0.0,
    )

    solve!(solver, (t0, tf); dt=dt, save_history=false)

    capTl = something(solver.cache.modelT.cap1_slab)
    capTs = something(solver.cache.modelT.cap2_slab)
    capCl = something(solver.cache.modelC.cap1_slab)
    capCs = something(solver.cache.modelC.cap2_slab)

    eTl = profile_error_on_phase(solver.state.Tlω, (x, t) -> case.Tl_exact(x, tf), capTl, tf)
    eTs = profile_error_on_phase(solver.state.Tsω, (x, t) -> case.Ts_exact(x, tf), capTs, tf)
    eCl = profile_error_on_phase(solver.state.Clω, (x, t) -> case.Cl_exact(x, tf), capCl, tf)
    eCs = profile_error_on_phase(solver.state.Csω, (x, t) -> case.Cs_exact(x, tf), capCs, tf)

    return (
        h=h,
        L2=(Tl=eTl.L2, Ts=eTs.L2, Cl=eCl.L2, Cs=eCs.L2),
        Linf=(Tl=eTl.Linf, Ts=eTs.Linf, Cl=eCl.Linf, Cs=eCs.Linf),
        totalL2=eTl.L2 + eTs.L2 + eCl.L2 + eCs.L2,
        solver=solver,
    )
end

@testset "MMS fixed-interface (1D)" begin
    o33 = _run_mms_fixed_1d(33)
    o49 = _run_mms_fixed_1d(49)
    o65 = _run_mms_fixed_1d(65)

    for out in (o33, o49, o65)
        @test isfinite(out.L2.Tl)
        @test isfinite(out.L2.Ts)
        @test isfinite(out.L2.Cl)
        @test isfinite(out.L2.Cs)
        @test isfinite(out.Linf.Tl)
        @test isfinite(out.Linf.Ts)
        @test isfinite(out.Linf.Cl)
        @test isfinite(out.Linf.Cs)
    end

    @test o49.totalL2 <= 1.02 * o33.totalL2
    @test o65.totalL2 <= 1.02 * o49.totalL2
    p12 = log(o33.totalL2 / o49.totalL2) / log(o33.h / o49.h)
    p23 = log(o49.totalL2 / o65.totalL2) / log(o49.h / o65.h)
    @test p12 > 0.9
    @test p23 > 0.9

    alg = interface_algebra_residuals(o65.solver)
    @test alg.continuity.maxabs < 1e-8
    @test alg.partition.maxabs < 1e-8
    @test alg.liquidus.maxabs < 1e-8
end

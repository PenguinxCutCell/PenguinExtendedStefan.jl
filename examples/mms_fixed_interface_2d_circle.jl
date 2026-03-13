using CartesianGrids
using PenguinBCs
using PenguinExtendedStefan

function run_level(n)
    case = manufactured_fixed_interface_2d_case(; center=(0.05, -0.03), radius=0.37)
    t0 = 0.1
    tf = 0.14
    dt0 = 0.003

    grid = CartesianGrid((-1.0, -1.0), (1.0, 1.0), (n, n))
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
        1e12,
        AlloyEquilibrium(case.k, case.Tm, case.m_liquidus);
        kappa_l=1.0,
        kappa_s=0.9,
        D_l=0.1,
        D_s=0.04,
        source_T_l=case.source_T_l,
        source_T_s=case.source_T_s,
        source_C_l=case.source_C_l,
        source_C_s=case.source_C_s,
    )
    opts = ExtendedStefanOptions(; scheme=:CN, reinit=false, coupling_max_iter=10)
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

    dt = dt0
    success = false
    last_err = nothing
    while dt >= 5e-4 && !success
        try
            solve!(solver, (t0, tf); dt=dt, save_history=false)
            success = all(isfinite, solver.state.Tlω) && all(isfinite, solver.state.Clω)
        catch err
            last_err = err
            success = false
        end
        success || (dt *= 0.5)
    end
    success || error("MMS 2D run failed for n=$n (last error: $(last_err))")

    capTl = something(solver.cache.modelT.cap1_slab)
    capTs = something(solver.cache.modelT.cap2_slab)
    capCl = something(solver.cache.modelC.cap1_slab)
    capCs = something(solver.cache.modelC.cap2_slab)

    eTl = profile_error_on_phase(solver.state.Tlω, (x, y, t) -> case.Tl_exact(x, y, tf), capTl, tf)
    eTs = profile_error_on_phase(solver.state.Tsω, (x, y, t) -> case.Ts_exact(x, y, tf), capTs, tf)
    eCl = profile_error_on_phase(solver.state.Clω, (x, y, t) -> case.Cl_exact(x, y, tf), capCl, tf)
    eCs = profile_error_on_phase(solver.state.Csω, (x, y, t) -> case.Cs_exact(x, y, tf), capCs, tf)

    return (n=n, dt=dt, totalL2=eTl.L2 + eTs.L2 + eCl.L2 + eCs.L2)
end

println("MMS fixed interface (2D circle)")
for n in (15, 23, 35)
    r = run_level(n)
    println("  n=$(r.n) dt=$(r.dt) total_L2=$(r.totalL2)")
end

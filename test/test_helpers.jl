using CartesianGrids
using PenguinBCs
using PenguinExtendedStefan

function basic_problem_1d(; n=65, s0=0.47, rhoL=6.0, kpart=0.7, mliq=-0.2, scheme=:BE)
    grid = CartesianGrid((0.0,), (1.0,), (n,))
    rep = LevelSetRep(grid, x -> x - s0)

    bcT = BorderConditions(; left=Dirichlet(1.0), right=Dirichlet(0.0))
    bcC = BorderConditions(; left=Dirichlet(1.0), right=Dirichlet(0.0))

    params = ExtendedStefanParams(
        1.0,
        1.0,
        rhoL,
        AlloyEquilibrium(kpart, 0.0, mliq);
        kappa_l=1.0,
        kappa_s=1.0,
        D_l=0.1,
        D_s=0.05,
        source_T_l=0.0,
        source_T_s=0.0,
        source_C_l=0.0,
        source_C_s=0.0,
    )

    opts = ExtendedStefanOptions(
        ;
        scheme=scheme,
        reinit=false,
        extend_iters=8,
        coupling_max_iter=10,
        coupling_tol=1e-8,
        coupling_reltol=1e-8,
        coupling_damping=0.7,
    )

    return ExtendedStefanProblem(grid, bcT, bcC, params, rep, opts)
end

function make_solver_from_problem(prob::ExtendedStefanProblem; t0=0.0)
    return build_solver(
        prob;
        t0=t0,
        Tlω0=0.4,
        Tlγ0=0.3,
        Tsω0=0.2,
        Tsγ0=0.3,
        Clω0=0.35,
        Clγ0=0.25,
        Csω0=0.2,
        Csγ0=0.175,
        speed0=0.0,
    )
end

function mk_reference_problem(n; t0=0.05, L=0.8)
    pars = martin_kauffman_reference()
    κT = pars.kappaT
    κS = pars.tau * κT
    h0 = martin_kauffman_interface_position(t0, pars)

    grid = CartesianGrid((-L,), (L,), (n,))
    rep = LevelSetRep(grid, x -> x - h0)

    bcT = BorderConditions(; left=Dirichlet(1.0), right=Dirichlet((x, t) -> martin_kauffman_fields(x, t, pars).T))
    bcC = BorderConditions(; left=Dirichlet(1.0), right=Dirichlet(0.0))

    params = ExtendedStefanParams(
        1.0,
        1.0,
        pars.S,
        AlloyEquilibrium(0.0, 0.0, -pars.Lambda);
        kappa_l=κT,
        kappa_s=κT,
        D_l=κS,
        D_s=1e-12,
        source_T_l=0.0,
        source_T_s=0.0,
        source_C_l=0.0,
        source_C_s=0.0,
    )

    opts = ExtendedStefanOptions(
        ;
        scheme=:BE,
        reinit=false,
        extend_iters=12,
        coupling_max_iter=16,
        coupling_tol=1e-8,
        coupling_reltol=1e-8,
        coupling_damping=0.7,
    )

    prob = ExtendedStefanProblem(grid, bcT, bcC, params, rep, opts)

    solver = build_solver(
        prob;
        t0=t0,
        Tlω0=(x, t) -> martin_kauffman_fields(x, t0, pars).T,
        Tlγ0=(x, t) -> martin_kauffman_fields(x, t0, pars).T,
        Tsω0=(x, t) -> martin_kauffman_fields(x, t0, pars).T,
        Tsγ0=(x, t) -> martin_kauffman_fields(x, t0, pars).T,
        Clω0=(x, t) -> martin_kauffman_fields(x, t0, pars).C,
        Clγ0=(x, t) -> martin_kauffman_fields(x, t0, pars).C,
        Csω0=0.0,
        Csγ0=0.0,
        speed0=pars.alpha * sqrt(κT / t0),
    )

    return solver, pars
end

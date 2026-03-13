using CartesianGrids
using PenguinBCs
using PenguinExtendedStefan

function run_level(n; t0=0.05, tf=0.08, L=0.8)
    pars = martin_kauffman_reference()
    κT = pars.kappaT

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
        D_l=pars.tau * κT,
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

    dx = 2L / (n - 1)
    dt = 0.2 * dx
    solve!(solver, (t0, tf); dt=dt, save_history=false)

    h_num = PenguinExtendedStefan.interface_position_from_levelset_1d(prob.interface_rep, grid)
    h_ref = martin_kauffman_interface_position(tf, pars)

    capTl = something(solver.cache.modelT.cap1_slab)
    capCl = something(solver.cache.modelC.cap1_slab)
    eT = profile_error_on_phase(solver.state.Tlω, (x, t) -> martin_kauffman_fields(x, tf, pars).T, capTl, tf).L2
    eC = profile_error_on_phase(solver.state.Clω, (x, t) -> martin_kauffman_fields(x, tf, pars).C, capCl, tf).L2

    return (n=n, dx=dx, e_h=abs(h_num - h_ref), eT=eT, eC=eC)
end

pars = martin_kauffman_reference()
println("1-D multicomponent melting (Martin and Kauffman, 1977)")
println("  S=$(pars.S), tau=$(pars.tau), Lambda=$(pars.Lambda)")
println("  A=$(pars.A), B=$(pars.B), alpha=$(pars.alpha)")

results = [run_level(n) for n in (17, 25, 33, 49, 65, 97)]
for r in results
    println("  n=$(r.n) dx=$(r.dx) e_h=$(r.e_h) eT=$(r.eT) eC=$(r.eC)")
end

ooc_h = [log(results[i].e_h / results[i+1].e_h) / log(results[i].dx / results[i+1].dx) for i in 1:length(results)-1]
ooc_T = [log(results[i].eT / results[i+1].eT) / log(results[i].dx / results[i+1].dx) for i in 1:length(results)-1]
ooc_C = [log(results[i].eC / results[i+1].eC) / log(results[i].dx / results[i+1].dx) for i in 1:length(results)-1]
println("  observed interface order: $(ooc_h)")
println("  observed temperature order: $(ooc_T)")
println("  observed concentration order: $(ooc_C)")

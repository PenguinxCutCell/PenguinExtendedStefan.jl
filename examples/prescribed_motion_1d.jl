using CartesianGrids
using PenguinBCs
using PenguinExtendedStefan

case = prescribed_planar_motion_case(; h0=0.36, V0=0.08, ω=0.0, A=0.0)

t0 = 0.0
tf = 0.06
dt = 0.003

grid = CartesianGrid((0.0,), (1.0,), (129,))
rep = LevelSetRep(grid, case.phi0)

bcT = BorderConditions(; left=Dirichlet((x, t) -> case.Tl_exact(x, t)), right=Dirichlet((x, t) -> case.Ts_exact(x, t)))
bcC = BorderConditions(; left=Dirichlet((x, t) -> case.Cl_exact(x, t)), right=Dirichlet((x, t) -> case.Cs_exact(x, t)))

params = ExtendedStefanParams(
    1.0,
    1.0,
    1.0,
    AlloyEquilibrium(0.7, 0.0, -0.2);
    kappa_l=1.0,
    kappa_s=0.8,
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

solve_prescribed_motion!(solver, (t0, tf); dt=dt, speed_callback=case.speed, save_history=false)

h_num = PenguinExtendedStefan.interface_position_from_levelset_1d(prob.interface_rep, prob.grid)
h_ref = case.h(tf)

println("Prescribed-motion MMS (1D)")
println("  interface position num = ", h_num)
println("  interface position ref = ", h_ref)
println("  |error|               = ", abs(h_num - h_ref))
println("  max |speed|           = ", maximum(abs.(vec(solver.state.speed_full))))

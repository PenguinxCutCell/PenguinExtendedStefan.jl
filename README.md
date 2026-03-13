# PenguinExtendedStefan.jl

`PenguinExtendedStefan.jl` is a sharp-interface binary-alloy Stefan solver in the PenguinxCutCell ecosystem.

It follows the `PenguinStefan.jl` workflow:
- level-set interface representation,
- `build_solver`, `step!`, `solve!` API,
- nonlinear fixed-point loop at each time step,
- one monolithic temperature+concentration linear solve per nonlinear iteration,
- interface speed update from the thermal Stefan closure.

## Governing Equations

Bulk equations (liquid/solid):

```math
\rho_l c_{p,l}\,\partial_t T_l = \nabla\cdot(\kappa_l\nabla T_l) + s_{T_l},
\qquad
\rho_s c_{p,s}\,\partial_t T_s = \nabla\cdot(\kappa_s\nabla T_s) + s_{T_s}
```

```math
\partial_t C_l = \nabla\cdot(D_l\nabla C_l) + s_{C_l},
\qquad
\partial_t C_s = \nabla\cdot(D_s\nabla C_s) + s_{C_s}
```

Normal convention: `n` points from solid to liquid.

Interface conditions:

```math
T_{l\Gamma}-T_{s\Gamma}=0
```
```math
C_{s\Gamma}-kC_{l\Gamma}=0
```
```math
T_{l\Gamma}-(T_m+mC_{l\Gamma})=0
```
```math
(D_s\partial_n C_s-D_l\partial_n C_l)-V_\Gamma(C_{l\Gamma}-C_{s\Gamma})=0
```

Speed closure (post-solve update):

```math
V_\Gamma = \frac{\kappa_s\partial_n T_s - \kappa_l\partial_n T_l}{\rho L}
```

## Relation To `PenguinStefan.jl`

`PenguinStefan.jl` targets classical thermal Stefan formulations.
`PenguinExtendedStefan.jl` extends that architecture to binary-alloy diffusion with partition/liquidus/solute-conservation coupling.

## Feature Matrix

| Feature | Status |
|---|---|
| Diphasic heat diffusion | implemented |
| Diphasic solute diffusion | implemented |
| Liquidus equilibrium | implemented |
| Partition law | implemented |
| Solute conservation | implemented |
| Stefan speed update | implemented |
| Level-set interface tracking | implemented |
| GHF | not implemented |
| Gibbs-Thomson | not implemented |
| Convection | not implemented |
| 3D validation | not implemented |

## Validation And Benchmark Suite

Validation ladder used in this package:

1. API/build and assembly smoke tests.
2. Fixed-interface interface-algebra checks (`TlΓ=TsΓ`, `CsΓ=kClΓ`, `TlΓ=Tm+mClΓ`).
3. Reduction-to-classical-Stefan checks (`k=1`, `m=0`, concentration decoupled).
4. 1D Martin-Kauffman benchmark (main binary-alloy reference).
5. Manufactured fixed-interface 1D and 2D cases.
6. Manufactured prescribed-motion 1D and 2D smoke checks.
7. Curved-interface 2D physics smoke checks.
8. Exploratory planar-front perturbation example (not a validated dendrite benchmark).

Primary benchmark reference:
- AFiD-MuRPhFi Stefan examples page, section
  [1-D multicomponent melting (Martin and Kauffman, 1977)](https://chowland.github.io/AFiD-MuRPhFi/examples/stefan/)

Exploratory instability context:
- Mullins-Sekerka planar-interface instability theory for dilute alloys.

## Quickstart

```julia
using CartesianGrids
using PenguinBCs
using PenguinExtendedStefan

grid = CartesianGrid((0.0,), (1.0,), (129,))
rep = LevelSetRep(grid, x -> x - 0.5)

bcT = BorderConditions(; left=Dirichlet(1.0), right=Dirichlet(0.0))
bcC = BorderConditions(; left=Dirichlet(1.0), right=Dirichlet(0.0))

params = ExtendedStefanParams(
    1.0,
    1.0,
    10.0,
    AlloyEquilibrium(0.7, 0.0, -0.2);
    kappa_l=1.0,
    kappa_s=1.0,
    D_l=0.1,
    D_s=0.05,
)

prob = ExtendedStefanProblem(grid, bcT, bcC, params, rep, ExtendedStefanOptions())
solver = build_solver(prob; Tlω0=0.5, Tlγ0=0.5, Tsω0=0.5, Tsγ0=0.5, Clω0=0.2, Clγ0=0.2, Csω0=0.1, Csγ0=0.1)
out = solve!(solver, (0.0, 0.1); dt=0.01)
metrics = extended_stefan_residual_metrics(out.solver)
```

## Examples

```bash
julia --project=. examples/fixed_interface_partition_liquidus_1d.jl
julia --project=. examples/reduction_to_classical_stefan_1d.jl
julia --project=. examples/martin_kauffman_1d_multicomponent_melting.jl
julia --project=. examples/mms_fixed_interface_1d.jl
julia --project=. examples/mms_fixed_interface_2d_circle.jl
julia --project=. examples/prescribed_motion_1d.jl
julia --project=. examples/prescribed_motion_2d_circle.jl
julia --project=. examples/planar_front_perturbation_2d.jl
julia --project=. examples/growth_instability_2d_animation.jl
```

## Current Limitations

Current version does not include Gibbs-Thomson, interface kinetics, global-height-function coupling, convection/Navier-Stokes, anisotropy, or phase-field formulations.

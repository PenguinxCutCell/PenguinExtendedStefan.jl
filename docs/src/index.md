# PenguinExtendedStefan.jl

`PenguinExtendedStefan.jl` is a sharp-interface binary-alloy Stefan package in the PenguinxCutCell ecosystem.

Implemented core features:
- diphasic heat diffusion,
- diphasic solute diffusion,
- liquidus equilibrium,
- partition law,
- solute conservation,
- thermal Stefan speed update,
- level-set interface tracking.

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
| Mushy zone | not implemented |
| 3D validation | not implemented |

## Workflow

The solver mirrors `PenguinStefan.jl`:
- `build_solver` creates the problem cache/state,
- `step!` advances one time step with nonlinear coupling iterations,
- `solve!` marches over a time span.

See:
- [Mathematical Model](model.md)
- [API](api.md)
- [Level-Set Workflow](levelset.md)
- [Validation](validation.md)
- [Examples](examples.md)
- [Diagnostics](diagnostics.md)
- [Limitations](limitations.md)

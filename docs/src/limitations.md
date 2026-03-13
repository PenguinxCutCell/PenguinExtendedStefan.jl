# Limitations / Next Steps

Current `v0.1` scope:
- level-set interface tracking only,
- diffusion-only heat/solute transport,
- monolithic `(T,C)` solve with lagged-speed solute row,
- thermal-Stefan speed update outside the linear solve.

Not implemented yet:
- global height function coupling,
- Gibbs-Thomson and kinetic interface laws,
- convection / Navier-Stokes coupling,
- mushy-region modeling,
- anisotropy,
- phase-field formulations,
- broad 3D validation.

Planned items are tracked in `TODO.md`.

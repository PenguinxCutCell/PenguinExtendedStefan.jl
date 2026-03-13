# Level-Set Workflow

`LevelSetRep` is the only interface representation in `v0.1`.

Per nonlinear iteration inside one time step:

1. Read `\phi^n` from the current interface state.
2. Use lagged speed field `V^{(m)}` to predict `\phi^{(m+1)}`.
3. Rebuild moving cut-cell slab geometry from `(\phi^n, \phi^{(m+1)})`.
4. Assemble and solve monolithic `T+C` system on that geometry.
5. Recompute interface speed from thermal Stefan flux balance.
6. Extend speed away from interface for robust level-set advection.
7. After convergence, advance `\phi` and optionally reinitialize.

Available hooks:
- `predict_phi`
- `extend_speed!`
- `advance!`
- `reinit!`

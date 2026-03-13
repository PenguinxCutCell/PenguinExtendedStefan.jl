# API

Primary exports:

- Interface representation:
  - `AbstractInterfaceRep`, `LevelSetRep`
- Parameter/options types:
  - `ExtendedStefanParams`, `ExtendedStefanOptions`
- Problem/state/solver:
  - `ExtendedStefanProblem`, `ExtendedStefanState`, `ExtendedStefanSolver`
- Solver entry points:
  - `build_solver`, `step!`, `solve!`
- Diagnostics:
  - `temperature_error_norms`
  - `concentration_error_norms`
  - `extended_stefan_residual_metrics`
  - `interface_algebra_residuals`
  - `stefan_residual_metrics`
  - `l2_error`
  - `linf_error`
  - `profile_error_on_phase`
  - `interface_position_error`
  - `phase_volume`, `active_mask`

- Validation helpers:
  - `martin_kauffman_reference`
  - `martin_kauffman_fields`
  - `martin_kauffman_interface_position`
  - `manufactured_fixed_interface_1d_case`
  - `manufactured_fixed_interface_2d_case`
  - `prescribed_planar_motion_case`
  - `prescribed_circle_motion_case`
  - `solve_prescribed_motion!`

Use REPL help for docstrings:

```julia
?ExtendedStefanProblem
?build_solver
?step!
```

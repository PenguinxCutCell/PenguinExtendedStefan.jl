# Diagnostics

`extended_stefan_residual_metrics(solver)` returns:
- thermal Stefan residual norm,
- solute conservation residual norm,
- speed increment norm,
- max speed,
- nonlinear iteration count.

`interface_algebra_residuals(solver)` returns interface-constraint residual norms for:
- temperature continuity,
- partition law,
- liquidus relation.

Field error helpers:
- `temperature_error_norms(...)`
- `concentration_error_norms(...)`
- `l2_error(...)`
- `linf_error(...)`
- `profile_error_on_phase(...)`
- `interface_position_error(...)`

Geometry reductions:
- `active_mask(cap)`
- `phase_volume(cap)`

These diagnostics are intended for regression tests, manufactured solutions, and mesh studies.

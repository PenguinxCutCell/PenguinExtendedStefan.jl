# Examples

Run from repository root:

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

Purpose by example:

- `fixed_interface_partition_liquidus_1d.jl`
  - checks interface algebra with nearly fixed interface.

- `reduction_to_classical_stefan_1d.jl`
  - demonstrates concentration-decoupled limit (`k=1`, `m=0`).

- `martin_kauffman_1d_multicomponent_melting.jl`
  - runs the Martin-Kauffman 1D multicomponent melting benchmark and prints interface-error convergence.

- `mms_fixed_interface_1d.jl`
  - manufactured fixed-interface verification in 1D.

- `mms_fixed_interface_2d_circle.jl`
  - manufactured fixed curved-interface verification in 2D.

- `prescribed_motion_1d.jl`
  - prescribed moving planar interface verification.

- `prescribed_motion_2d_circle.jl`
  - prescribed moving circle smoke verification.

- `planar_front_perturbation_2d.jl`
  - exploratory planar-front perturbation run for morphological-instability-onset studies.
  - this is not a validated dendrite benchmark.

- `growth_instability_2d_animation.jl`
  - writes a CairoMakie MP4 animation with side-by-side temperature and concentration fields.

using Test
using CartesianGrids
using PenguinBCs
using PenguinExtendedStefan

include("test_helpers.jl")

include("test_bc_types.jl")
include("test_layout_and_build.jl")
include("test_interface_algebra_1d.jl")
include("test_reduction_classical_stefan_1d.jl")
include("test_martin_kauffman_1d.jl")
include("test_mms_fixed_interface_1d.jl")
include("test_mms_fixed_interface_2d.jl")
include("test_prescribed_motion_1d.jl")
include("test_prescribed_motion_2d_smoke.jl")
include("test_curved_interface_2d_smoke.jl")

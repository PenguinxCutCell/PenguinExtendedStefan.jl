module PenguinExtendedStefan

using LinearAlgebra
using SparseArrays
using StaticArrays
using LinearSolve

using CartesianGrids
using CartesianGeometry
using CartesianOperators
using PenguinBCs
using PenguinSolverCore
using PenguinDiffusion
using LevelSetMethods
import PenguinSolverCore: solve!

export AbstractInterfaceRep
export LevelSetRep
export ExtendedStefanParams, ExtendedStefanOptions
export ExtendedStefanProblem, ExtendedStefanState, ExtendedStefanSolver
export build_solver, step!, solve!
export temperature_error_norms, concentration_error_norms, extended_stefan_residual_metrics
export phase_volume, active_mask
export l2_error, linf_error, profile_error_on_phase, interface_position_error
export interface_algebra_residuals, stefan_residual_metrics
export martin_kauffman_reference, martin_kauffman_fields, martin_kauffman_interface_position
export manufactured_fixed_interface_1d_case, manufactured_fixed_interface_2d_case
export prescribed_planar_motion_case, prescribed_circle_motion_case, solve_prescribed_motion!

include("interface/abstract.jl")
include("types.jl")
include("validation/reference_fields.jl")
include("validation/interface_algebra.jl")
include("validation/martin_kauffman_1977.jl")
include("validation/manufactured_fixed_interface.jl")
include("validation/manufactured_prescribed_motion.jl")
include("validation/manufactured_1d.jl")
include("validation/reductions.jl")
include("interface/levelset_rep.jl")
include("interface/predictors.jl")
include("moving/slab_geometry.jl")
include("assembly/cache.jl")
include("assembly/interface_rows.jl")
include("assembly/monolithic_system.jl")
include("speed/stefan_speed.jl")
include("stepper.jl")
include("diagnostics.jl")
include("io.jl")

end # module

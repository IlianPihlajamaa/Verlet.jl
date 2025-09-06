module Constraints

using LinearAlgebra
using ..Core
include("constraints.jl")

# Export DistanceConstraints type
export DistanceConstraints, apply_rattle!, apply_shake!, velocity_verlet_shake_rattle!, constraint_residuals, remove_com_motion!

end

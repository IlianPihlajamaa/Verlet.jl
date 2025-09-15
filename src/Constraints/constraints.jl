module Constraints

using ..Core
using LinearAlgebra, StaticArrays

include("constraints.jl")

export DistanceConstraints, apply_shake!, apply_rattle!, velocity_verlet_shake_rattle!, remove_com_motion!, constraint_residuals

end

module Verlet
using LinearAlgebra, StaticArrays, StaticArrays

# Import submodules
include("Core/Core.jl")
include("Potentials/Potentials.jl")
include("Neighbors/Neighbors.jl")
include("Constraints/Constraints.jl")
include("Thermostats/Thermostats.jl")

using .Core
using .Neighbors
using .Potentials
using .Constraints
using .Thermostats

# Re-export public API for backward compatibility
export ParticleSystem, velocity_verlet!, kinetic_energy, potential_energy
export CubicBox, minimum_image!, wrap_positions!, box_length
export lj_forces
export NeighborList, build_neighborlist
export DistanceConstraints, apply_shake!
export langevin_baoab!, degrees_of_freedom, instantaneous_temperature, velocity_rescale!

export maybe_rebuild!, build_cellgrid, build_neighborlist_cells, rebin!
export DistanceConstraints, apply_rattle!, apply_shake!, velocity_verlet_shake_rattle!, constraint_residuals, remove_com_motion!
export langevin_baoab_constrained!
end # module Verlet

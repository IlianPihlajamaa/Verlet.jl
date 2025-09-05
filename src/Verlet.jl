module Verlet
using LinearAlgebra

# Public API exports (implementations live in separate files)
export ParticleSystem, velocity_verlet!, kinetic_energy, potential_energy
export langevin_baoab_constrained!, constraint_residuals
export CubicBox, minimum_image!, wrap_positions!, lj_forces
export NeighborList, build_neighborlist, maybe_rebuild!, max_displacement_since_build
export CellGrid, build_cellgrid, rebin!, build_neighborlist_cells
export langevin_baoab!, instantaneous_temperature, degrees_of_freedom, velocity_rescale!
export DistanceConstraints, velocity_verlet_shake_rattle!, apply_shake!, apply_rattle!, remove_com_motion!

# Implementation files
include("particles.jl")
include("integrators.jl")
include("boxes.jl")
include("forces.jl")
include("cellgrid.jl")
include("neighborlist_cells.jl")
include("neighborlist.jl")
include("constraints.jl")
include("thermostats.jl")

end # module Verlet

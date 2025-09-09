"""
Core module: fundamental types & utilities.

Data representation (current):
  positions :: Vector{SVector{D,T}}
  velocities :: Vector{SVector{D,T}}
  forces(R) :: Vector{SVector{D,T}} (or (F,U) with return_potential=true)

Matrix-based NxD layouts are deprecated.
"""
module Core
export T_float, T_int, Dims

# Default type selectors
default_float() = Float64
default_int() = Int
const T_float = default_float()
const T_int = default_int()
const Dims = 3  # spatial dimensionality (default, can be overridden)

export ParticleSystem, velocity_verlet!, kinetic_energy, potential_energy
export CubicBox, minimum_image, wrap_positions!, box_length

include("particles.jl")
include("integrators.jl")
include("forces.jl")
include("boxes.jl")

end

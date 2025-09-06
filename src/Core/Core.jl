module Core
export T_float, T_int, Dims

# Default type selectors
default_float() = Float64
default_int() = Int
const T_float = default_float()
const T_int = default_int()
const Dims = 3  # spatial dimensionality (default, can be overridden)
const Dims = 3  # spatial dimensionality

export ParticleSystem, velocity_verlet!, kinetic_energy, potential_energy
export CubicBox, minimum_image!, wrap_positions!, box_length

include("particles.jl")
include("integrators.jl")
include("forces.jl")
include("boxes.jl")
include("utils.jl")

end

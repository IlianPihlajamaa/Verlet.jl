module Core

using LinearAlgebra, StaticArrays

const T_Float = Float64
const T_Int = Int64
const Dims = 3

include("boxes.jl")
include("particles.jl")
include("integrator_interface.jl")
include("energy.jl")
include("forces.jl")
include("potential_types.jl")
include("neighbor_types.jl")

export System, natoms, natomtypes, AbstractBox, AbstractIntegrator, integrate!, potential_energy, kinetic_energy
export CubicBox, minimum_image, wrap_positions!, box_length
export T_Float, T_Int, Dims
export PotentialNeighborList
export AbstractPotentialPair, AbstractPairPotential, AbstractBondPotential, AbstractAnglePotential, AbstractDihedralPotential, AbstractImproperPotential, AbstractNeighborList
export compute_forces!, compute_all_forces!, ForceField

end

module Core

using LinearAlgebra, StaticArrays

const T_Float = Float64
const T_int = Int64
const Dims = 3

include("boxes.jl")
include("particles.jl")
include("integrators.jl")
include("forces.jl")
include("potential_types.jl")
include("neighbor_types.jl")

export System, natoms, natomtypes, AbstractBox, velocity_verlet!, potential_energy, kinetic_energy
export CubicBox, minimum_image, wrap_positions!, box_length
export T_Float, T_int, Dims
export NeighborPair, PotentialNeighborList, MasterNeighborEntry, MasterNeighborList
export AbstractPotentialPair, AbstractPairPotential, AbstractBondPotential, AbstractAnglePotential, AbstractDihedralPotential, AbstractImproperPotential

end

module Core

using LinearAlgebra, StaticArrays, StructArrays

const T_Float = Float64
const T_int = Int64
const Dims = 3

include("potentials.jl")
include("boxes.jl")
include("particles.jl")
include("integrators.jl")
include("forces.jl")
include("neighbor_types.jl")
include("bonded_types.jl")
include("bonded_forces.jl")

export System, natoms, natomtypes, AbstractBox, velocity_verlet!, potential_energy, kinetic_energy, compute_forces!
export AbstractBondPotential, HarmonicBond, Bond, AbstractAnglePotential, HarmonicAngle, Angle, AbstractDihedralPotential, PeriodicDihedral, Dihedral
export CubicBox, minimum_image, wrap_positions!, box_length
export T_Float, T_int, Dims
export AbstractPotentialPair, AbstractPairPotential, AbstractBondPotential
export NeighborPair, PotentialNeighborList, MasterNeighborEntry, MasterNeighborList

end

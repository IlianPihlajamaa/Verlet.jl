module Potentials

using ..Core

# Potentials - defines different potential energy functions

# Abstract types and structs for pairwise potentials
include("potentials.jl")

# Specific pairwise potentials
include("lj.jl")
include("coulomb.jl")

# Bonded potentials
include("bonded.jl")

export LennardJones, Coulomb, LJPair, CoulPair, PairTable, Bond, Angle, Dihedral, HarmonicBond, HarmonicAngle, PeriodicDihedral, compute_forces!, potential_energy
export AbstractPotentialPair, AbstractPairPotential, AbstractBondPotential, AbstractAnglePotential, AbstractDihedralPotential, AbstractImproperPotential
export NeighborPair, PotentialNeighborList, MasterNeighborEntry, MasterNeighborList

end

module Potentials

import ..Core
using ..Neighbors

# Potentials - defines different potential energy functions

# Abstract types and structs for pairwise potentials
include("pairpotentials.jl")

# Specific pairwise potentials
include("lj.jl")
include("coulomb.jl")

# Bonded potentials
include("bonded.jl")

export LennardJones, Coulomb, LJPair, CoulPair, PairTable, Bond, Angle, Dihedral, HarmonicBond, HarmonicAngle, PeriodicDihedral,  potential_energy
export  AbstractBondPotential, AbstractAnglePotential, AbstractDihedralPotential, AbstractImproperPotential
export NeighborPair, PotentialNeighborList, MasterNeighborEntry, MasterNeighborList

end

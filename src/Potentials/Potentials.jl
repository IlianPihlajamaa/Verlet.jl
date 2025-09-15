module Potentials
using LinearAlgebra, StaticArrays, StructArrays

using ..Core
using ..Neighbors

# Potentials - defines different potential energy functions

# Abstract types and structs for pairwise potentials
include("pair_potentials.jl")

# Specific pairwise potentials
include("lj.jl")
include("coulomb.jl")

# Bonded potentials
include("bonded.jl")

export LennardJones, Coulomb, LJPair, CoulPair, PairTable, Bond, Angle, Dihedral, HarmonicBond, HarmonicAngle, PeriodicDihedral, compute_forces!, potential_energy

end

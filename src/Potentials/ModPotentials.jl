module Potentials
using LinearAlgebra, StaticArrays, StructArrays

using ..Core
using ..Neighbors

include("pairpotentials.jl")
include("lj.jl")
include("coulomb.jl")

# Bonded potentials
include("bonded.jl")

export LennardJones, Coulomb, LJPair, CoulPair, PairTable, Bond, Angle, Dihedral, HarmonicBond, HarmonicAngle, PeriodicDihedral, compute_forces!, potential_energy

end
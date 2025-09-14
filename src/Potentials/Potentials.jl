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

export LJPair, CoulPair, PairTable, LennardJones, Coulomb, Bond, Angle, Dihedral, HarmonicBond, HarmonicAngle, PeriodicDihedral, compute_forces!

end

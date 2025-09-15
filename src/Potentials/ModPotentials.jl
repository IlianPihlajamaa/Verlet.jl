module Potentials
using LinearAlgebra, StaticArrays, StructArrays

using ..Core
using ..Neighbors

<<<<<<< HEAD
include("pairpotentials.jl")
=======
# Potentials - defines different potential energy functions

# Abstract types and structs for pairwise potentials
include("pairpotentials.jl")

# Specific pairwise potentials
>>>>>>> e020923 (fix all namespaces)
include("lj.jl")
include("coulomb.jl")

# Bonded potentials
include("bonded.jl")

export LennardJones, Coulomb, LJPair, CoulPair, PairTable, Bond, Angle, Dihedral, HarmonicBond, HarmonicAngle, PeriodicDihedral, compute_forces!, potential_energy

<<<<<<< HEAD
end
=======
end
>>>>>>> e020923 (fix all namespaces)

abstract type AbstractPotentialPair end
abstract type AbstractPairPotential end

# Using StructArrays for efficient storage and access
using StructArrays

"""
    LJPair{T<:Number}

Parameters for a Lennard-Jones interaction between two particle types.
"""
struct LJPair{T<:Number} <: AbstractPotentialPair
    ε::T
    σ::T
    rc::T
end

"""
    CoulPair{T<:Number}

Parameters for a Coulomb interaction between two particle types.
"""
struct CoulPair{T<:Number} <: AbstractPotentialPair
    q1q2::T # product of charges
    rc::T
end

"""
    PairTable{F<:AbstractPotentialPair}

A table to store parameters for each pair of particle types.
"""
struct PairTable{F<:AbstractPotentialPair}
    table::Matrix{F}
end

"""
    NeighborPair{F<:AbstractPotentialPair, IntT<:Integer}

A pair of interacting particles `i` and `j` with their interaction parameters.
"""
struct NeighborPair{F<:AbstractPotentialPair, IntT<:Integer}
    i::IntT
    j::IntT
    pair::F
end

"""
    PotentialNeighborList{F<:AbstractPotentialPair}

A neighbor list for a specific potential, storing pairs and their parameters.
It is a `StructArray` of `NeighborPair`s.
"""
const PotentialNeighborList{F} = StructArray{NeighborPair{F, T_int}} where {F<:AbstractPotentialPair}


"""
    MasterNeighborEntry

An entry in the master neighbor list, containing the indices of a pair of particles
and the squared distance between them at the time of construction.
"""
struct MasterNeighborEntry
    i::T_int
    j::T_int
    r2::T_Float
end

"""
    MasterNeighborList{T<:Number}

The master neighbor list, containing all candidate pairs within the largest cutoff radius.
"""
struct MasterNeighborList{T<:Number}
    skin::T
    entries::Vector{MasterNeighborEntry}
    # Other fields like reference positions could be added here if needed
end

using ..Core
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


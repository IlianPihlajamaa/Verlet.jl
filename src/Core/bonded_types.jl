# Abstract types for bonded potentials
abstract type AbstractBondPotential end
abstract type AbstractAnglePotential end
abstract type AbstractDihedralPotential end
abstract type AbstractImproperPotential end

# Concrete implementations of bonded potentials
struct HarmonicBond{T} <: AbstractBondPotential
    k::T    # Spring constant
    r0::T   # Equilibrium distance
end

struct HarmonicAngle{T} <: AbstractAnglePotential
    k::T    # Spring constant
    θ0::T   # Equilibrium angle
end

struct PeriodicDihedral{T} <: AbstractDihedralPotential
    k::T    # Barrier height
    n::Int  # Periodicity
    ϕ0::T   # Phase offset
end

# Structs to hold the particle indices for each interaction
# These are what will be stored in the System's specific_potentials list
struct Bond
    i::Int
    j::Int
    potential::AbstractBondPotential
end

struct Angle
    i::Int
    j::Int
    k::Int
    potential::AbstractAnglePotential
end

struct Dihedral
    i::Int
    j::Int
    k::Int
    l::Int
    potential::AbstractDihedralPotential
end

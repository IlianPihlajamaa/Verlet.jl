
# Implementation for particle state, e.g. ParticleSystem, CubicBox
using StaticArrays, LinearAlgebra



"""
    ParticleSystem{Dims,T_float}

A lightweight container for particle states.

# Fields
- `positions::Vector{SVector{Dims,T_float}}`: length-N vector of positions (each an SVector)
- `velocities::Vector{SVector{Dims,T_float}}`: length-N vector of velocities
- `masses::Vector{T_float}`: length-N vector of particle masses

"""
mutable struct ParticleSystem{Dims,T_float}
    positions::Vector{SVector{Dims,T_float}}   # (N)
    velocities::Vector{SVector{Dims,T_float}}  # (N)
    masses::Vector{T_float}                    # (N)
end

"""
    kinetic_energy(system::ParticleSystem) -> T_float

Total kinetic energy: `∑ ½ mᵢ ‖vᵢ‖²`.
"""
function kinetic_energy(system::ParticleSystem{Dims,T_float}) where {Dims,T_float}
    @assert length(system.positions) == length(system.velocities) "positions/velocities must be same size"
    @assert length(system.masses) == length(system.positions) "length(masses) must equal number of particles"
    Ekin = zero(T_float)

    for i in 1:length(system.positions)
        Ekin += T_float(0.5) * system.masses[i] * dot(system.velocities[i], system.velocities[i])
    end
    return Ekin
end

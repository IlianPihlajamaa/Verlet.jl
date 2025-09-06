# Implementation for particle state, e.g. ParticleSystem, CubicBox


"""
    ParticleSystem{T}

A lightweight container for particle states.

# Fields
- `positions::Matrix{T}`: `(N × D)` array of positions (rows = particles).
- `velocities::Matrix{T}`: `(N × D)` array of velocities.
- `masses::Vector{T}`: length-`N` vector of particle masses.

!!! tip "Memory layout"
    For best performance, keep `positions` and `velocities` as `Matrix{T}` with fixed `D` (1–3 typical).
"""
mutable struct ParticleSystem{T_float}
    positions::Matrix{T_float}   # (N × D)
    velocities::Matrix{T_float}  # (N × D)
    masses::Vector{T_float}      # (N)
end

"""
    kinetic_energy(system::ParticleSystem) -> T_float

Total kinetic energy: `∑ ½ mᵢ ‖vᵢ‖²`.
"""
function kinetic_energy(system::ParticleSystem{T_float}) where {T_float}
    @assert size(system.positions) == size(system.velocities) "positions/velocities must be same size"
    @assert length(system.masses) == size(system.positions, 1) "length(masses) must equal number of particles"
    v2 = sum(abs2, system.velocities; dims=2)
    return T_float(0.5) * sum(system.masses .* vec(v2))
end

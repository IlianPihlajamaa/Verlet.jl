# Implementation for particle state, e.g. ParticleSystem, CubicBox

"""
    ParticleSystem

A lightweight container for particle states.

# Fields
- `positions::Matrix{Float64}`: `(N × D)` array of positions (rows = particles).
- `velocities::Matrix{Float64}`: `(N × D)` array of velocities.
- `masses::Vector{Float64}`: length-`N` vector of particle masses.

!!! tip "Memory layout"
    For best performance, keep `positions` and `velocities` as `Matrix{Float64}` with fixed `D` (1–3 typical).
"""
mutable struct ParticleSystem
    positions::Matrix{Float64}   # (N × D)
    velocities::Matrix{Float64}  # (N × D)
    masses::Vector{Float64}      # (N)
end

"""
    kinetic_energy(system::ParticleSystem) -> Float64

Total kinetic energy: `∑ ½ mᵢ ‖vᵢ‖²`.
"""
function kinetic_energy(system::ParticleSystem)::Float64
    @assert size(system.positions) == size(system.velocities) "positions/velocities must be same size"
    @assert length(system.masses) == size(system.positions, 1) "length(masses) must equal number of particles"
    v2 = sum(abs2, system.velocities; dims=2)
    return 0.5 * sum(system.masses .* vec(v2))
end

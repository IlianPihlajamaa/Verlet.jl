
"""
    module Verlet

Basic velocity Verlet integrator for tiny molecular-dynamics style problems.

# Overview

- **Data layout:** Positions and velocities are `N×D` `Matrix{Float64}` (rows = particles, columns = spatial dimensions).
- **Masses:** Length-`N` `Vector{Float64}`.
- **Forces:** User-supplied function that takes positions `(N×D)` and returns forces `(N×D)` in the same units.

# Units

You pick units, but they must be **consistent**:
- Positions `r` in length (e.g., m, Å),
- Velocities `v` in length/time,
- Masses `m` in mass,
- Forces `F` in mass·length/time²,
- Timestep `dt` in time.

# Stability & Pitfalls

- Choose `dt` small enough relative to the fastest vibrational frequency in your system.
- Use `Float64` throughout to minimize drift.
- Energy is approximately conserved for small `dt`, but will drift if `dt` is too large.

# Exports

[`ParticleSystem`](@ref), [`velocity_verlet!`](@ref), [`kinetic_energy`](@ref), [`potential_energy`](@ref)
"""
module Verlet
export ParticleSystem, velocity_verlet!, kinetic_energy, potential_energy

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
    velocity_verlet!(system::ParticleSystem, forces::Function, dt::Float64)

Advance `system` by one timestep `dt` using the velocity Verlet integrator.

The `forces` function **must** accept a positions matrix `(N × D)` and return an
`(N × D)` matrix of forces (same shape).

!!! info "Algorithm"
    Given positions `r`, velocities `v`, accelerations `a = F(r)/m`, and `dt`:
    1. `r ← r + v*dt + 0.5*a*dt^2`
    2. Recompute `a' = F(r)/m`
    3. `v ← v + 0.5*(a + a')*dt`

# Examples
```julia
julia> using Verlet

julia> forces(r) = zeros(size(r))  # free particle
forces (generic function with 1 method)

julia> ps = ParticleSystem([0.0 0.0], [1.0 0.0], [1.0]);

julia> velocity_verlet!(ps, forces, 0.1);  # advances in-place

julia> ps.positions[1,1]
0.1
````

!!! warning "Shape rules"
- `size(positions) == size(velocities)`
- `length(masses) == size(positions, 1)`
"""
 function velocity_verlet!(system::ParticleSystem, forces::Function, dt::Float64)
    @assert size(system.positions) == size(system.velocities) "positions/velocities must be same size"
    @assert length(system.masses) == size(system.positions, 1) "length(masses) must equal number of particles"
    dt2 = dt * dt
     # Current accelerations a = F/m, broadcasting over rows (per particle)
     a = forces(system.positions) ./ system.masses
     # Update positions
        system.positions .+= system.velocities .* dt .+ 0.5 .* a .* dt2
     # New accelerations from updated positions
    a_new = forces(system.positions) ./ system.masses
    # Update velocities
    system.velocities .+= 0.5 .* (a .+ a_new) .* dt
    return system
end

"""
    kinetic_energy(system::ParticleSystem) -> Float64

Total kinetic energy: `∑ ½ mᵢ ‖vᵢ‖²`.

# Examples

```julia
julia> ps = ParticleSystem([0.0 0.0; 1.0 0.0],
                           [1.0 2.0; 0.0 0.0],
                           [1.0, 2.0]);

julia> kinetic_energy(ps)
0.5*(1*(1^2+2^2) + 2*(0^2+0^2))  # = 2.5
2.5
```

"""
function kinetic_energy(system::ParticleSystem)::Float64
    @assert size(system.positions) == size(system.velocities) "positions/velocities must be same size"
    @assert length(system.masses) == size(system.positions, 1) "length(masses) must equal number of particles"
    # Sum velocities^2 per particle (row); result is N×1, so vec to N-vector
    v2 = sum(abs2, system.velocities; dims=2)
    return 0.5 * sum(system.masses .* vec(v2))
end

"""
    potential_energy(system::ParticleSystem, forces::Function) -> Float64

Try to obtain **total potential energy** from the user-supplied `forces` function.

Supported convention:

* If `forces(r; return_potential=true)` is supported, it must return `(F, U)`
  where `F` is the `(N × D)` force matrix and `U` is a scalar total potential energy.

If this convention is not supported, an error is thrown.

# Examples

```julia
julia> function ho_forces(r; return_potential=false)
           F = -r
           return return_potential ? (F, 0.5 * sum(abs2, r)) : F
       end;

julia> ps = ParticleSystem([1.0 0.0], [0.0 0.0], [1.0]);

julia> potential_energy(ps, ho_forces)
0.5
```

"""
function potential_energy(system::ParticleSystem, forces::Function)::Float64
    # Try convention 2: keyword to request potential
    F_U_kw = try
        forces(system.positions; return_potential=true)
    catch
        nothing
    end
    if F_U_kw !== nothing && F_U_kw isa Tuple && length(F_U_kw) == 2
        _, U = F_U_kw
        return float(U)
    end
   error("Force function does not support `return_potential=true`; cannot compute potential_energy.")
end


end # module Verlet

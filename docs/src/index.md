# Verlet.jl

A minimal **velocity Verlet** integrator for tiny MD-style problems.

!!! note
    You choose the units; just keep them consistent across positions `r`, velocities `v`, masses `m`,
    forces `F`, and timestep `dt`.

## Quickstart

```@example quickstart
using Verlet

# Free particle in 2D
forces(r) = zeros(size(r))

ps = ParticleSystem([0.0 0.0],  # 1×2 positions
                    [1.0 0.0],  # 1×2 velocities
                    [1.0])      # masses

dt = 0.1
velocity_verlet!(ps, forces, dt)
ps.positions
```

## Harmonic oscillator
```@example ho
using Verlet

# Hooke's law with k = 1, potential U = 0.5 * |r|^2
function ho_forces(r; return_potential=false)
    F = -r
    return return_potential ? (F, 0.5 * sum(abs2, r)) : F
end

ps = ParticleSystem([1.0 0.0], [0.0 0.0], [1.0])

dt = 0.1
for _ in 1:100
    velocity_verlet!(ps, ho_forces, dt)
end

(kin = kinetic_energy(ps), pot = potential_energy(ps, ho_forces))
```

## Energy monitoring

```@example energy
using Verlet

function ho_forces(r; return_potential=false)
    F = -r
    return return_potential ? (F, 0.5 * sum(abs2, r)) : F
end

ps = ParticleSystem([1.0 0.0], [0.0 1.0], [1.0])
dt = 0.05

energies = Float64[]
for _ in 1:200
    velocity_verlet!(ps, ho_forces, dt)
    push!(energies, kinetic_energy(ps) + potential_energy(ps, ho_forces))
end

(round(minimum(energies), digits=6), round(maximum(energies), digits=6))
```

## Performance tips

* Keep arrays as `Matrix{Float64}` / `Vector{Float64}` to avoid type instability.
* Prefer **in-place** force computations in your own code paths; if you must allocate, reuse buffers.
* Avoid huge `dt`. Start small (e.g., `1e-3` in your time units) and increase cautiously.

See also: \[Numerics & Pitfalls]\(@ref numerics).
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

## Speeding up LJ with a Neighbor List

The naive Lennard–Jones force kernel scales as **O(N²)**, which quickly becomes expensive as the number of particles grows.
A **Verlet neighbor list** reduces this to **O(N)** on average at fixed density by storing nearby pairs and updating them only occasionally.

```@example
using Verlet

# Create a cubic periodic box
box = CubicBox(20.0)

# Random positions for 100 particles in 3D
R = randn(100, 3)
wrap_positions!(R, box)

# Build neighbor list with cutoff + skin
nlist = build_neighborlist(R, box; cutoff=2.5, skin=0.4)

# Rebuild when particles have moved far enough
maybe_rebuild!(nlist, R, box) && @info "Neighbor list rebuilt"

# Compute forces using the list-aware kernel
F, U = lj_forces(R, box, nlist; rcut=2.5, return_potential=true)
@show size(F), U
```

### API

```@docs
NeighborList
build_neighborlist
maybe_rebuild!
max_displacement_since_build
wrap_positions!
```

### Performance Tips

- Choose a **skin** large enough that rebuilds are infrequent, but not so large that each particle has too many neighbors.  
    A good starting point is `skin ≈ 0.3σ`.
- Ensure the box length `L` satisfies `L > 2*(cutoff + skin)` to avoid ambiguous minimum-image distances.
- Accumulated forces and energies are stored in `Float64` for stability even if input positions are `Float32`.
- Rebuilds are **O(N²)**, but are triggered rarely. Per-step force computation with the neighbor list is **O(N)** on average.

### Common Pitfalls

- Using a **too small skin** can cause missed interactions if particles move across the buffer before a rebuild.
- Forgetting to call [`wrap_positions!`](@ref) regularly may lead to large apparent displacements across periodic boundaries.
- The neighbor list includes symmetric neighbors. Forces should only be applied once per pair (the built-in `lj_forces` handles this).

### See Also

- [`build_neighborlist`](@ref) — construct a new list from positions.  
- [`maybe_rebuild!`](@ref) — keep a list up-to-date based on displacements.  
- [`wrap_positions!`](@ref) — enforce periodic wrapping of coordinates.  
- [`lj_forces`](@ref) — compute Lennard–Jones forces with or without a neighbor list.

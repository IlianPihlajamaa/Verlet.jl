# Verlet.jl

A minimal **velocity Verlet** integrator for tiny MD-style problems.

!!! note
    You choose the units; just keep them consistent across positions `r`, velocities `v`, masses `m`,
    forces `F`, and timestep `dt`.

## Quickstart

```@example quickstart
using Verlet, StaticArrays

# Free particle in 2D
positions = [@SVector [0.0, 0.0]]
velocities = [@SVector [1.0, 0.0]]
masses = [1.0]
forces(R) = [@SVector zeros(2) for _ in R]

ps = ParticleSystem(positions, velocities, masses)

dt = 0.1
velocity_verlet!(ps, forces, dt)
ps.positions
```

## Next Steps

Check out the \[Guide → Constrained Dynamics](@ref constraints-guide) section to learn how to:

- Set up bond constraints with [`DistanceConstraints`](@ref)
- Run constrained dynamics with [`velocity_verlet_shake_rattle!`](@ref)


## Harmonic oscillator

```@example ho
using Verlet, StaticArrays, LinearAlgebra

# Hooke's law with k = 1, potential U = 0.5 * |r|^2
function ho_forces(R; return_potential=false)
    F = [-r for r in R]
    U = 0.5 * sum(norm(r)^2 for r in R)
    return return_potential ? (F, U) : F
end

positions = [@SVector [1.0, 0.0]]
velocities = [@SVector [0.0, 0.0]]
masses = [1.0]
ps = ParticleSystem(positions, velocities, masses)

dt = 0.1
for _ in 1:100
    velocity_verlet!(ps, ho_forces, dt)
end

(kin = kinetic_energy(ps), pot = potential_energy(ps, ho_forces))
```

## Energy monitoring


```@example energy
using Verlet, StaticArrays, LinearAlgebra

function ho_forces(R; return_potential=false)
    F = [-r for r in R]
    U = 0.5 * sum(norm(r)^2 for r in R)
    return return_potential ? (F, U) : F
end

positions = [@SVector [1.0, 0.0]]
velocities = [@SVector [0.0, 1.0]]
masses = [1.0]
ps = ParticleSystem(positions, velocities, masses)
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

With constraints, you can simulate rigid bonds (e.g. water models) and safely
increase timestep sizes while preserving stability.

## Speeding up LJ with a Neighbor List

The naive Lennard–Jones force kernel scales as **O(N²)**, which quickly becomes expensive as the number of particles grows.
A **Verlet neighbor list** reduces this to **O(N)** on average at fixed density by storing nearby pairs and updating them only occasionally.

```@example
using Verlet

# Create a cubic periodic box
box = CubicBox(20.0)

using StaticArrays
# Random positions for 100 particles in 3D as SVectors
R = [@SVector randn(3) for _ in 1:100]
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

See the full API reference in `api.md`.

## Note on Particle Representation

Verlet.jl now uses `Vector{SVector{Dims, T_Float}}` (from StaticArrays) to represent particle positions, velocities, displacements, and forces. This provides better performance and type stability compared to the previous `Matrix`-based approach. All user code and force functions should now expect and return vectors of SVectors, e.g.:

```@example 1
using Verlet, StaticArrays
N, D = 100, 3
positions = [@SVector randn(D) for _ in 1:N]
velocities = [@SVector zeros(D) for _ in 1:N]
masses = ones(N)
ps = ParticleSystem(positions, velocities, masses)
```

Force functions should accept and return `Vector{SVector}` as well:

```julia
function my_forces(R)
    # R is Vector{SVector{D, T}}
    return [@SVector zeros(length(R[1])) for _ in R]
end
```

## NEW: O(N) Build with Cell-Linked Lists + Half Neighbor Lists

The classic `build_neighborlist` uses an **O(N²)** construction. For larger systems
you can switch to a **cell-linked grid** builder that is **O(N)** at fixed density
and emits a **half list** (each pair stored once with `j > i`):

```@example 1
D = 3
box = CubicBox(10.0)
cutoff, skin = 2.5, 0.4
R = [SVector{D}((rand(D) .- 0.5) .* box.L) for _ in 1:2_000]  # random positions in (-L/2, L/2]
wrap_positions!(R, box)
grid = build_cellgrid(R, box; cell_size=cutoff+skin)
nl = build_neighborlist_cells(R, box; cutoff=cutoff, skin=skin, grid=grid)
# Force evaluation with half list (branch-free inner loop)
F = lj_forces(R, box, nl; rcut=cutoff)  # or (F,U) with return_potential=true
size(F)
```

### Why half lists?

\* **Less memory**: store each pair once (≈½ the entries of a symmetric list).
\* **Fewer branches**: kernel no longer checks `j > i` at runtime.
\* **Same physics**: forces are accumulated for both `i` and `j` when a pair is visited.

### Rebuild policy

Use the same **half-skin rule** via `maybe_rebuild!`. With O(N) builds you
can afford a **smaller skin** (e.g., 0.2–0.3) to reduce neighbor count and tighten
force errors:

```@example 1
box = CubicBox(15.0)
R = randn(SVector{3, Float64}, 500); wrap_positions!(R, box)
cutoff, skin = 2.5, 0.3
grid = build_cellgrid(R, box; cell_size=cutoff+skin)
nl = build_neighborlist_cells(R, box; cutoff=cutoff, skin=skin, grid=grid)

# ... advance dynamics updating R ...

if maybe_rebuild!(nl, R, box)
    rebin!(grid, R, box) # O(N)
    nl = build_neighborlist_cells(R, box; cutoff=cutoff, skin=skin, grid=grid)
end
```

### Pitfalls

\* **Box too small**: ensure `L > 2*(cutoff + skin)` so minimum-image distances are unambiguous.
\* **Cell size semantics**: the grid uses an **effective** width `L/nx ≥ cutoff+skin`.
\* **Units**: `R`/`L` must share the same units; `cutoff`/`skin` are in those units.
\* **Precision**: distance math uses `Float64`. Mixed-precision inputs are converted.

### Performance Tips

- Choose a **skin** large enough that rebuilds are infrequent, but not so large that each particle has too many neighbors.  
    A good starting point is `skin ≈ 0.3σ`.
- Ensure the box length `L` satisfies `L > 2*(cutoff + skin)` to avoid ambiguous minimum-image distances.
- Accumulated forces and energies are stored in `Float64` for stability even if input positions are `Float32`.
- Rebuilds are **O(N²)**, but are triggered rarely. Per-step force computation with the neighbor list is **O(N)** on average.

### Common Pitfalls

- Using a **too small skin** can cause missed interactions if particles move across the buffer before a rebuild.
- Forgetting to call `wrap_positions!` regularly may lead to large apparent displacements across periodic boundaries.
- The neighbor list includes symmetric neighbors. Forces should only be applied once per pair (the built-in `lj_forces` handles this).

### See Also

- [`build_neighborlist`](@ref) — construct a new list from positions.  
- [`maybe_rebuild!`](@ref) — keep a list up-to-date based on displacements.  
- [`wrap_positions!`](@ref) — enforce periodic wrapping of coordinates.  
- [`lj_forces`](@ref) — compute Lennard–Jones forces with or without a neighbor list.

## Further Notes

* [Constraints](@ref) — how SHAKE/RATTLE are applied and how to monitor residuals
* [Numerical Notes](@ref) — guidance on tolerances, thermostat interaction, and reproducibility

---

```@contents
Pages = ["constraints.md", "numerics.md", "api.md"]
Depth = 2
```

```@index
Pages = ["api.md"]
```

# Verlet.jl

A minimal **velocity Verlet** integrator for tiny MD-style problems.

!!! note
    You choose the units; just keep them consistent across positions `r`, velocities `v`, masses `m`,
    forces `F`, and timestep `dt`.

## Quickstart

```@example quickstart
using Verlet, StaticArrays

# 1. Set up the system
positions = [@SVector [0.0, 0.0, 0.0]]
velocities = [@SVector [1.0, 0.0, 0.0]]
forces_storage = [@SVector [0.0, 0.0, 0.0]]
masses = [1.0]
box = CubicBox(10.0)
types = [1]
type_names = Dict(1 => :A)
sys = System(positions, velocities, forces_storage, masses, box, types, type_names)

# 2. Define a potential (e.g., Lennard-Jones)
# Note: For a single particle, the force will be zero. This is just for demonstration.
ϵ = 1.0
σ = 1.0
rc = 2.5
lj_pair = Verlet.Potentials.LJPair(ϵ, σ, rc)
params = Verlet.Potentials.PairTable(fill(lj_pair, (1, 1)))
exclusions = Tuple{Verlet.Core.T_int,Verlet.Core.T_int}[]
lj = Verlet.Potentials.LennardJones(params, exclusions, 0.5)
ff = Verlet.Neighbors.ForceField((lj,))

# 3. Define a force function compatible with the integrator
function compute_forces_for_integrator(positions, system, forcefield, master_nl)
    system.positions .= positions # Update positions in the system object
    Verlet.Neighbors.build_all_neighbors!(master_nl, forcefield, system)
    Verlet.Neighbors.compute_all_forces!(system, forcefield)
    return system.forces
end

# 4. Run the simulation
dt = 0.1
master_nl = Verlet.Neighbors.MasterNeighborList(0.5)
# Wrap the force function to match the integrator's signature
force_wrapper(R) = compute_forces_for_integrator(R, sys, ff, master_nl)
Verlet.Core.velocity_verlet!(sys, force_wrapper, dt)
sys.positions
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

positions = [@SVector [1.0, 0.0, 0.0]]
velocities = [@SVector [0.0, 0.0, 0.0]]
forces = [@SVector [0.0, 0.0, 0.0]]
masses = [1.0]
box = CubicBox(10.0)
types = [1]
type_names = Dict(1 => :A)
sys = System(positions, velocities, forces, masses, box, types, type_names)

dt = 0.1
for _ in 1:100
    velocity_verlet!(sys, ho_forces, dt)
end

(pot = ho_forces(sys.positions, return_potential=true)[2])
```

## Energy monitoring


```@example energy
using Verlet, StaticArrays, LinearAlgebra

function ho_forces(R; return_potential=false)
    F = [-r for r in R]
    U = 0.5 * sum(norm(r)^2 for r in R)
    return return_potential ? (F, U) : F
end

positions = [@SVector [1.0, 0.0, 0.0]]
velocities = [@SVector [0.0, 1.0, 0.0]]
forces = [@SVector [0.0, 0.0, 0.0]]
masses = [1.0]
box = CubicBox(10.0)
types = [1]
type_names = Dict(1 => :A)
sys = System(positions, velocities, forces, masses, box, types, type_names)
dt = 0.05

energies = Float64[]
for _ in 1:200
    velocity_verlet!(sys, ho_forces, dt)
    push!(energies, ho_forces(sys.positions, return_potential=true)[2]) # kinetic_energy removed
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

## ForceField API for Potentials

The recommended way to handle pair potentials like Lennard-Jones is with the `ForceField` API. This provides a flexible way to combine multiple potentials and uses an efficient neighbor list implementation.

```@example
using Verlet, StaticArrays

# 1. Set up a system
box = CubicBox(20.0)
R = [@SVector randn(3) for _ in 1:100]
wrap_positions!(R, box)
sys = System(
    R,
    [@SVector(zeros(3)) for _ in R],
    [@SVector(zeros(3)) for _ in R],
    ones(length(R)),
    box,
    ones(Int, length(R)),
    Dict(1 => :A)
)

# 2. Define a Lennard-Jones potential
lj = Verlet.Potentials.LennardJones(
    Verlet.Potentials.PairTable(fill(Verlet.Potentials.LJPair(1.0, 1.0, 2.5), (1, 1))),
    Tuple{Verlet.Core.T_int,Verlet.Core.T_int}[],
    0.5
)

# 3. Create a ForceField
ff = Verlet.Neighbors.ForceField((lj,))

# 4. Build neighbor lists and compute forces
master_skin = 0.5
master_nl = Verlet.Neighbors.MasterNeighborList(master_skin)
Verlet.Neighbors.build_all_neighbors!(master_nl, ff, sys)
Verlet.Neighbors.compute_all_forces!(sys, ff)

@show sys.forces[1]
```

### Neighbor List Methods

You can choose the neighbor list algorithm with the `method` keyword in `build_all_neighbors!`:
- `:cells` (default): Fast `O(N)` cell-based algorithm.
- `:bruteforce`: Slower `O(N^2)` algorithm for debugging.
- `:all_pairs`: Includes all pairs, ignoring cutoffs.

```julia
master_nl = Verlet.Neighbors.MasterNeighborList(master_skin)
Verlet.Neighbors.build_all_neighbors!(master_nl, ff, sys, method=:bruteforce)
```

- [`build_master_neighborlist!`](@ref) — construct a new master neighbor list.
- [`wrap_positions!`](@ref) — enforce periodic wrapping of coordinates.

# 4. Run the simulation
dt = 0.1
# Wrap the force function to match the integrator's signature
force_wrapper(R) = compute_forces_for_integrator(R, sys, ff)
velocity_verlet!(sys, force_wrapper, dt)
sys.positions
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

positions = [@SVector [1.0, 0.0, 0.0]]
velocities = [@SVector [0.0, 0.0, 0.0]]
forces = [@SVector [0.0, 0.0, 0.0]]
masses = [1.0]
box = CubicBox(10.0)
types = [1]
type_names = Dict(1 => :A)
sys = System(positions, velocities, forces, masses, box, types, type_names)

dt = 0.1
for _ in 1:100
    velocity_verlet!(sys, ho_forces, dt)
end

(pot = ho_forces(sys.positions, return_potential=true)[2])
```

## Energy monitoring


```@example energy
using Verlet, StaticArrays, LinearAlgebra

function ho_forces(R; return_potential=false)
    F = [-r for r in R]
    U = 0.5 * sum(norm(r)^2 for r in R)
    return return_potential ? (F, U) : F
end

positions = [@SVector [1.0, 0.0, 0.0]]
velocities = [@SVector [0.0, 1.0, 0.0]]
forces = [@SVector [0.0, 0.0, 0.0]]
masses = [1.0]
box = CubicBox(10.0)
types = [1]
type_names = Dict(1 => :A)
sys = System(positions, velocities, forces, masses, box, types, type_names)
dt = 0.05

energies = Float64[]
for _ in 1:200
    velocity_verlet!(sys, ho_forces, dt)
    push!(energies, ho_forces(sys.positions, return_potential=true)[2]) # kinetic_energy removed
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

## ForceField API for Potentials

The recommended way to handle pair potentials like Lennard-Jones is with the `ForceField` API. This provides a flexible way to combine multiple potentials and uses an efficient neighbor list implementation.

```@example
using Verlet, StaticArrays

# 1. Set up a system
box = CubicBox(20.0)
R = [@SVector randn(3) for _ in 1:100]
wrap_positions!(R, box)
sys = System(
    R,
    [@SVector(zeros(3)) for _ in R],
    [@SVector(zeros(3)) for _ in R],
    ones(length(R)),
    box,
    ones(Int, length(R)),
    Dict(1 => :A)
)

# 2. Define a Lennard-Jones potential
lj = LennardJones(
    PairTable(fill(LJPair(1.0, 1.0, 2.5), (1, 1))),
    Tuple{T_int,T_int}[],
    0.5
)

# 3. Create a ForceField
ff = ForceField((lj,))

# 4. Build neighbor lists and compute forces
master_skin = 0.5
build_all_neighbors!(ff, sys, master_skin)
compute_all_forces!(sys, ff)

@show sys.forces[1]
```

### Neighbor List Methods

You can choose the neighbor list algorithm with the `method` keyword in `build_all_neighbors!`:
- `:cells` (default): Fast `O(N)` cell-based algorithm.
- `:bruteforce`: Slower `O(N^2)` algorithm for debugging.
- `:all_pairs`: Includes all pairs, ignoring cutoffs.

```julia
build_all_neighbors!(ff, sys, master_skin, method=:bruteforce)
```

- [`build_master_neighborlist`](@ref) — construct a new master neighbor list.
- [`wrap_positions!`](@ref) — enforce periodic wrapping of coordinates.

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

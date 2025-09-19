# Tutorial 2 · Pair Potentials & Neighbour Lists

This tutorial expands the quickstart example into a small Lennard-Jones fluid.
You will learn how to parameterise pair potentials, let `ForceField` manage
neighbour lists, and monitor potential energy during a simulation.

## Goal

Simulate 64 identical particles interacting via a Lennard-Jones (12-6) pair
potential inside a cubic periodic box.

## 1. Assemble a particle system

```@example lj
using Random, StaticArrays, Verlet

Random.seed!(2024)
N = 64
box = CubicBox(8.0)
randv() = SVector{3}(randn(), randn(), randn())
positions = [randv() for _ in 1:N]
wrap_positions!(positions, box)
zero3 = SVector{3}(0.0, 0.0, 0.0)
velocities = fill(zero3, N)
forces = fill(zero3, N)
masses = fill(1.0, N)
types = ones(Int, N)
type_names = Dict(1 => :A)
```

## 2. Build a Lennard-Jones potential

We store `ε`, `σ`, and cutoff `r_c` in a `PairTable`. The per-type pair object
is looked up automatically using the particle type IDs.

```@example lj
ϵ  = 1.0
σ  = 1.0
rc = 2.5

params = PairTable(fill(LJPair(ϵ, σ, rc), (1, 1)))
exclusions = Tuple{T_Int, T_Int}[]
lj = LennardJones(params, exclusions, 0.5)
ff = ForceField((lj,))
```

## 3. Construct the system with forces attached

```@example lj
sys = System(positions, velocities, forces, masses, box, types, type_names;
             forcefield = ff)
```

The `Neighbors.ForceField` stores an internal master neighbour list. It lazily
initialises that list the first time `compute_all_forces!` (or `integrate!`) is
invoked.

## 4. Run dynamics and monitor energy

We run unconstrained velocity Verlet for a short trajectory and capture the
potential energy after each step.

```@example lj
vv = VelocityVerlet(5e-3)
energies = Float64[]
integrate!(vv, sys, 200; callback = (system, step, _) -> begin
    push!(energies, compute_potential_energy(system))
    return nothing
end)

minimum(energies), maximum(energies)
```

`compute_potential_energy(system)` automatically rebuilds neighbour lists when
particles move more than half the `skin` distance since the last build.

## 5. Inspect neighbours explicitly

You can trigger a manual neighbour rebuild and examine the stored pairs:

```@example lj
master = MasterNeighborList(sys; cutoff = rc, skin = lj.skin)
Verlet.Neighbors.build_all_neighbors!(master, ff, sys)
length(master.pairs)
```

Each pair is an `SVector{2,Int}` containing `(i, j)` indices. Concrete pair
potentials maintain their own neighbour lists derived from this master object.

## Recap & next steps

- Pair potentials are parameterised via `PairTable`s and added to a
  `Neighbors.ForceField`.
- Neighbour lists rebuild automatically during `compute_all_forces!` and can be
  inspected or forced manually.
- Energies can be sampled in a callback during `integrate!`.

Continue with [Tutorial 3 · Constraints in Practice](constraints.md) to enforce
rigid bonds during the dynamics.

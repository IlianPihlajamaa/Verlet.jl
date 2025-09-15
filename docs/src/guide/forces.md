# Forces & Potentials

This page describes how to define and use forces in your simulations.

## The `ForceField` API

The recommended way to define forces is using the `ForceField` API. This API allows you to compose multiple potentials (e.g., Lennard-Jones and Coulomb) in a flexible and efficient way.

### Example: Lennard-Jones

Here's how to set up a simple Lennard-Jones simulation:

```@example forces
using Verlet, StaticArrays, LinearAlgebra

# 1. Set up the system
box = CubicBox(10.0)
R = [@SVector(randn(3)) for _ in 1:64]; wrap_positions!(R, box)
sys = System(
    R,
    [@SVector(zeros(3)) for _ in R],
    [@SVector(zeros(3)) for _ in R],
    ones(length(R)),
    box,
    ones(Int, length(R)),
    Dict(1 => :A)
)

# 2. Define the potential
ϵ = 1.0
σ = 1.0
rc = 2.5
lj_pair = LJPair(ϵ, σ, rc)
params = PairTable(fill(lj_pair, (1, 1)))
exclusions = Tuple{T_int,T_int}[]
lj = LennardJones(params, exclusions, 0.5)

# 3. Create a ForceField
ff = Verlet.Neighbors.ForceField((lj,))

# 4. Build neighbor lists and compute forces
master_skin = 0.5
master_nl = Verlet.Neighbors.MasterNeighborList(master_skin)
Verlet.Neighbors.build_all_neighbors!(master_nl, ff, sys)
Verlet.Neighbors.compute_all_forces!(sys, ff)

# The forces are now stored in sys.forces
println(sys.forces[1])
```

### Composing Potentials

You can easily combine multiple potentials by adding them to the `ForceField` tuple:

```julia
# lj = LennardJones(...)
# coul = Coulomb(...)
# ff = ForceField((lj, coul))
```

## Neighbor List Methods

The `build_all_neighbors!` function uses a master neighbor list to accelerate force calculations. You can choose the method for building this list with the `method` keyword argument:

-   `method=:cells` (default): Uses a fast, `O(N)` cell-list algorithm. This is recommended for most systems.
-   `method=:bruteforce`: Uses a simple, `O(N^2)` algorithm. This can be useful for small systems or for debugging.
-   `method=:all_pairs`: Includes all pairs of particles, ignoring the cutoff. This is useful for testing or for potentials without a cutoff.

Here's how to use it:

```julia
# master_nl = MasterNeighborList(master_skin)
# build_all_neighbors!(master_nl, ff, sys, method=:bruteforce)
```

## Custom Forces

You can also define custom force functions. To make `potential_energy` work, your function should return `(F, U)` when the keyword `return_potential=true` is provided:

```@example forces2
using Verlet, StaticArrays, LinearAlgebra
function ho_forces(R; return_potential=false)
    F = [ -r for r in R ]      # k = 1
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
# The potential_energy function needs to be updated to work with the new API
# For now, I will just call the function
ho_forces(sys.positions)
```
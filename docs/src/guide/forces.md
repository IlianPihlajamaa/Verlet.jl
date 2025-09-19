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
exclusions = Tuple{T_Int,T_Int}[]
lj = LennardJones(params, exclusions, 0.5)

# 3. Create a ForceField
ff = Verlet.Neighbors.ForceField((lj,))

# 4. Compute forces (neighbor lists are managed by the ForceField)
Verlet.Neighbors.compute_all_forces!(sys, ff)

# The forces are now stored in sys.forces
println(sys.forces[1])
```

### Composing Potentials

You can easily combine multiple potentials by adding them to the `ForceField` tuple:

```julia
lj = LennardJones(...)
coul = Coulomb(...)
ff = ForceField((lj, coul))
```

## Neighbor List Methods

The `build_all_neighbors!` function uses a master neighbor list to accelerate force calculations. You can choose the method for building this list with the `method` keyword argument:

-   `method=:cells` (default): Uses a fast, `O(N)` cell-list algorithm. This is recommended for most systems.
-   `method=:bruteforce`: Uses a simple, `O(N^2)` algorithm. This can be useful for small systems or for debugging.
-   `method=:all_pairs`: Includes all pairs of particles, ignoring the cutoff. This is useful for testing or for potentials without a cutoff.

Here's how to use it:

```julia
master_nl = MasterNeighborList(sys; cutoff=rc, skin=master_skin)
build_all_neighbors!(master_nl, ff, sys, method=:bruteforce)
```

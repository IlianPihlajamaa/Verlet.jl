# Verlet.jl

A minimal **velocity Verlet** integrator for tiny MD-style problems.

!!! note
    You choose the units; just keep them consistent across positions `r`, velocities `v`, masses `m`,
    forces `F`, and timestep `dt`.

## Quickstart

```@example quickstart
using Verlet, StaticArrays

struct Hooke
    k::Float64
end

function Verlet.Core.compute_forces!(pot::Hooke, sys::System)
    @inbounds for i in 1:natoms(sys)
        sys.forces[i] += -pot.k * sys.positions[i]
    end
    return sys
end

box = CubicBox(10.0)
positions = [@SVector [1.0, 0.0, 0.0]]
velocities = [@SVector [0.0, 0.0, 0.0]]
forces_storage = [@SVector [0.0, 0.0, 0.0]]
masses = [1.0]
types = [1]
type_names = Dict(1 => :A)
ff = ForceField((Hooke(1.0),))
sys = System(positions, velocities, forces_storage, masses, box, types, type_names; forcefield=ff)

vv = VelocityVerlet(0.05)
integrate!(vv, sys, 100)
sys.positions
```



## Harmonic oscillator

```@example ho
using Verlet, StaticArrays, LinearAlgebra

struct HO
    k::Float64
end

function Verlet.Core.compute_forces!(pot::HO, sys::System)
    @inbounds for i in 1:natoms(sys)
        sys.forces[i] += -pot.k * sys.positions[i]
    end
    return sys
end

function ho_energy(sys::System, pot::HO)
    0.5 * pot.k * sum(norm(r)^2 for r in sys.positions)
end

pot = HO(1.0)
positions = [@SVector [1.0, 0.0, 0.0]]
velocities = [@SVector [0.0, 0.0, 0.0]]
forces = [@SVector [0.0, 0.0, 0.0]]
masses = [1.0]
box = CubicBox(10.0)
types = [1]
type_names = Dict(1 => :A)
sys = System(positions, velocities, forces, masses, box, types, type_names;
           forcefield=ForceField((pot,)))

vv = VelocityVerlet(0.1)
integrate!(vv, sys, 100)

ho_energy(sys, pot)
```

## Energy monitoring


```@example energy
using Verlet, StaticArrays, LinearAlgebra

struct HO
    k::Float64
end

function Verlet.Core.compute_forces!(pot::HO, sys::System)
    @inbounds for i in 1:natoms(sys)
        sys.forces[i] += -pot.k * sys.positions[i]
    end
    return sys
end

ho_energy(sys::System, pot::HO) = 0.5 * pot.k * sum(norm(r)^2 for r in sys.positions)

pot = HO(1.0)
positions = [@SVector [1.0, 0.0, 0.0]]
velocities = [@SVector [0.0, 1.0, 0.0]]
forces = [@SVector [0.0, 0.0, 0.0]]
masses = [1.0]
box = CubicBox(10.0)
types = [1]
type_names = Dict(1 => :A)
sys = System(positions, velocities, forces, masses, box, types, type_names;
           forcefield=ForceField((pot,)))
dt = 0.05

energies = Float64[]
vv = VelocityVerlet(dt)
integrate!(vv, sys, 200;
           callback = (sys, step, _) -> begin
               push!(energies, kinetic_energy(sys) + ho_energy(sys, pot))
               return nothing
           end)

(round(minimum(energies), digits=6), round(maximum(energies), digits=6))
```


## ForceField API for Potentials

The recommended way to handle pair potentials like Lennard-Jones is with the `ForceField` API. This provides a flexible way to combine multiple potentials and uses an efficient neighbor list implementation.

```@example
using Verlet, StaticArrays

# 1. Set up a system
box = CubicBox(20.0)
R = [@SVector randn(3) for _ in 1:100]
wrap_positions!(R, box)

# 2. Define a Lennard-Jones potential
rc = 2.5
lj = Verlet.Potentials.LennardJones(
    Verlet.Potentials.PairTable(fill(Verlet.Potentials.LJPair(1.0, 1.0, rc), (1, 1))),
    Tuple{Verlet.Core.T_Int,Verlet.Core.T_Int}[],
    0.5
)

# 3. Create a ForceField and system
ff = Verlet.Neighbors.ForceField((lj,))
sys = System(
    R,
    [@SVector(zeros(3)) for _ in R],
    [@SVector(zeros(3)) for _ in R],
    ones(length(R)),
    box,
    ones(Int, length(R)),
    Dict(1 => :A);
    forcefield=ff
)

# 4. Compute forces — the ForceField maintains its master neighbor list
Verlet.Neighbors.compute_all_forces!(sys, ff)

@show sys.forces[1]

dt = 0.01
vv = VelocityVerlet(dt)
integrate!(vv, sys, 10)
```

### Neighbor List Methods

You can choose the neighbor list algorithm with the `method` keyword in `build_all_neighbors!`:
- `:cells` (default): Fast `O(N)` cell-based algorithm.
- `:bruteforce`: Slower `O(N^2)` algorithm for debugging.
- `:all_pairs`: Includes all pairs, ignoring cutoffs.

```julia
rc = 2.5
master_nl = Verlet.Neighbors.MasterNeighborList(sys; cutoff=rc, skin=master_skin)
Verlet.Neighbors.build_all_neighbors!(master_nl, ff, sys, method=:bruteforce)
```

- [`build_master_neighborlist!`](@ref) — construct a new master neighbor list.
- [`wrap_positions!`](@ref) — enforce periodic wrapping of coordinates.

## Next Steps

Check out the \[Guide → Constrained Dynamics](@ref constraints-guide) section to learn how to:

- Set up bond constraints with [`DistanceConstraints`](@ref)
- Run constrained dynamics with [`velocity_verlet_shake_rattle!`](@ref)


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

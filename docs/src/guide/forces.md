# Forces & Potentials

This guide explains how forces plug into Verlet.jl, from bespoke
`compute_forces!` implementations to the built-in Lennard-Jones and Coulomb
potentials managed by `Neighbors.ForceField`.

## 1. Bring your own potential

Any concrete type becomes a potential once you implement
`Verlet.Core.compute_forces!(pot, system)`. The method must accumulate forces
into `system.forces` and return the system.

```@example forces
using Random, StaticArrays, Verlet

Random.seed!(3)

struct Springs
    k::Float64
end

function Verlet.Core.compute_forces!(pot::Springs, sys::System)
    @inbounds for i in 1:natoms(sys)
        sys.forces[i] += -pot.k * sys.positions[i]
    end
    return sys
end

function Verlet.Core.compute_potential_energy(pot::Springs, sys::System)
    energy = zero(eltype(sys.positions[1]))
    @inbounds for i in 1:natoms(sys)
        energy += 0.5 * pot.k * sum(abs2, sys.positions[i])
    end
    return energy
end
```

Add the potential to a `ForceField` so it participates in the standard force
pipeline:

```@example forces
springs = Springs(0.25)
ff = ForceField((springs,))
```

During `compute_all_forces!(sys, ff)` the force accumulator is zeroed, all layers
in `ff.layers` are evaluated, followed by any `sys.specific_potentials`.

## 2. Pair potentials with neighbour lists

Verlet.jl ships Lennard-Jones and Coulomb interactions that look up parameters
from a `PairTable` based on particle types. The `Neighbors.ForceField` wrapper
manages master neighbour lists behind the scenes.

```@example forces
N = 32
box = CubicBox(6.0)
randv() = SVector{3}(randn(), randn(), randn())
positions = [randv() for _ in 1:N]
wrap_positions!(positions, box)
zero3 = SVector{3}(0.0, 0.0, 0.0)
velocities = fill(zero3, N)
forces = fill(zero3, N)
masses = fill(1.0, N)
types = ones(Int, N)
type_names = Dict(1 => :A)

params = PairTable(fill(LJPair(1.0, 1.0, 2.5), (1, 1)))
lj = LennardJones(params, Tuple{T_Int,T_Int}[], 0.4)
ff = ForceField((springs, lj))

sys = System(positions, velocities, forces, masses, box, types, type_names;
             forcefield = ff)
Verlet.Core.compute_all_forces!(sys)  # neighbour lists prepared automatically
sys.forces[1]
```

Combine multiple pair potentials simply by stacking them inside the tuple. Each
layer gets its own per-potential neighbour list derived from the shared master
list.

## 3. Bonded interactions

Bonds, angles, and dihedrals operate on explicit particle indices and live in
the `specific_potentials` tuple on the `System`.

```@example forces
bond = Bond(1, 2, HarmonicBond(100.0, 1.0))
angle = Angle(1, 2, 3, HarmonicAngle(50.0, deg2rad(109.5)))

data = System(positions, velocities, forces, masses, box, types, type_names;
              forcefield = ff, specific_potentials = (bond, angle))
Verlet.Core.compute_all_forces!(data)
```

Bonded terms run after every layer listed in the `ForceField`, ensuring their
contributions accumulate alongside pair forces.

## 4. Energy evaluation

The same building blocks compute potential energy. Supply the forcefield via
`system.forcefield` or pass it explicitly:

```@example forces
E_pairs = compute_potential_energy(sys)           # uses sys.forcefield
E_pairs == compute_potential_energy(sys, ff)
```

If you provide a standalone force callback that supports
`return_potential=true`, the helper `Verlet.Core.potential_energy(system,
forces)` extracts the second element of the returned tuple.

## 5. When to rebuild neighbours manually

Most workflows rely on the automatic rebuild triggered by `compute_all_forces!`.
You only need manual control when modifying positions without calling the
integrator (e.g. Monte Carlo moves). In that case:

```@example forces
master = MasterNeighborList(sys; cutoff = 2.5, skin = 0.4)
Verlet.Neighbors.build_all_neighbors!(master, ff, sys)
```

Follow up with `maybe_rebuild(sys, master)` whenever particle displacements may
exceed `skin/2`.

For a tutorial-style walkthrough, see
[Tutorial 2 · Pair Potentials & Neighbour Lists](../tutorials/pair_potentials.md).

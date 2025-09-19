# Building Systems

A `Verlet.Core.System` bundles together particle data (positions, velocities,
forces, masses, types) and the objects required to compute interactions. This
page summarises the key fields, constructors, and utility helpers.

## Core fields

The `System` type stores positions, velocities, forces, masses, particle types,
and optional forcefield/bonded interaction hooks. See the API reference for the
full constructor signature.

Each vector of particle data must have the same length. The default element
types use `Float64` and `Int`, but any concrete floating-point and integer types
can be employed.

!!! tip "StaticVectors recommended"
    Positions, velocities, and forces are stored as `Vector{SVector{D,T}}`. This
    choice keeps per-particle operations allocation-free and type-stable. The
    `StaticArrays.jl` package offers convenient constructors.

## Constructing a system

```@example system
using Random, StaticArrays, Verlet

Random.seed!(1)
box = CubicBox(10.0)
randv() = SVector{3}(randn(), randn(), randn())
positions  = [randv() for _ in 1:4]
zero3 = SVector{3}(0.0, 0.0, 0.0)
velocities = fill(zero3, 4)
forces     = fill(zero3, 4)
masses     = fill(1.0, 4)
types      = ones(Int, 4)
type_names = Dict(1 => :A)

sys = System(positions, velocities, forces, masses, box, types, type_names)
```

No forces are computed yet—the `forcefield` slot defaults to `nothing`. Attach a
force field to take advantage of automated neighbour management:

```@example system
lj = LennardJones(PairTable(fill(LJPair(1.0, 1.0, 2.5), (1, 1))), Tuple{T_Int,T_Int}[], 0.4)
ff = ForceField((lj,))
sys = System(positions, velocities, forces, masses, box, types, type_names;
             forcefield = ff)
```

## Particle types & metadata

- `types::Vector{IT}` stores integer labels per particle. Use them to look up
  pair parameters from `PairTable`s.
- `type_names::Dict{IT,Symbol}` maps labels to human-friendly symbols (handy for
  output or assigning new parameters).
- Add your own metadata by storing side tables keyed by particle index—`System`
does not attempt to own every attribute.

## Bonded interactions

To add bonds, angles, or dihedrals, populate `specific_potentials` with the
corresponding interaction objects. They are evaluated after the force-field
layers during `compute_all_forces!` and `compute_potential_energy`.

```@example system
bond = Bond(1, 2, HarmonicBond(100.0, 1.0))
angle = Angle(1, 2, 3, HarmonicAngle(50.0, deg2rad(109.5)))

sys = System(positions, velocities, forces, masses, box, types, type_names;
             specific_potentials = (bond, angle))
```

## Utilities

- `natoms(sys)` and `natomtypes(sys)` report particle counts.
- `kinetic_energy(sys)` returns the total kinetic energy.
- `wrap_positions!(R, box)` keeps coordinates inside the primary simulation cell.
- `maybe_rebuild(system, master::MasterNeighborList)` checks whether a master
   neighbour list needs rebuilding (call before force evaluations if managing
   neighbours yourself).

## Common patterns

- **Allocate once:** reuse the same `System` when scanning parameters or
  performing temperature ramps—mutate the existing vectors instead of creating
  new ones each iteration.
- **Custom fields:** if you require additional arrays (e.g. per-particle
  dipoles), maintain them separately and keep indices aligned with the core
  vectors.

Continue with [Forces & Potentials](forces.md) to learn how custom interactions
hook into the forcefield machinery.

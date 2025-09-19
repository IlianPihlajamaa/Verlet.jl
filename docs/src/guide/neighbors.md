# Neighbour Lists

Efficient pair-force evaluation relies on culling distant particles. Verlet.jl
bundles a master neighbour list together with per-potential lists populated on
demand. This page summarises the main types and workflows.

## Master neighbour list

`MasterNeighborList` stores candidate pairs within `cutoff + skin`, caching buffers so rebuilds stay inexpensive.
Create one from either a `System` or raw positions plus box. The list stores all
candidate pairs within `cutoff + skin` and reuses buffers across rebuilds.

```@example neighbors
using Random, StaticArrays, Verlet

Random.seed!(2)
N = 16
box = CubicBox(6.0)
randv() = SVector{3}(randn(), randn(), randn())
R = [randv() for _ in 1:N]
wrap_positions!(R, box)
master = MasterNeighborList(R, box; cutoff = 2.5, skin = 0.4)
length(master.pairs)
```

## ForceField integration

`Neighbors.ForceField` extends `Verlet.Core.prepare_neighbors!` so that
`compute_all_forces!` rebuilds the master list only when particles move more
than `skin/2` from their reference positions.

```@example neighbors
params = PairTable(fill(LJPair(1.0, 1.0, 2.5), (1, 1)))
lj = LennardJones(params, Tuple{T_Int,T_Int}[], 0.4)
ff = ForceField((lj,))
```

When you call `compute_all_forces!(sys, ff)` or integrate using a `System` whose
`forcefield` is `ff`, neighbour lists are prepared automatically.

## Manual rebuilds

Sometimes you may want explicit control over neighbour maintenance, for example
in Monte Carlo moves. Use `build_all_neighbors!` with a pre-allocated master
list:

```@example neighbors
zero3 = SVector{3}(0.0, 0.0, 0.0)
sys = System(R, fill(zero3, N), fill(zero3, N),
             fill(1.0, N), box, ones(Int, N), Dict(1 => :A); forcefield = ff)

Verlet.Neighbors.build_all_neighbors!(master, ff, sys)
length(master.pairs)
```

`maybe_rebuild(sys, master)` checks the maximum squared displacement since the
last rebuild and triggers a rebuild only if needed:

```@example neighbors
master_before = length(master.pairs)
# Perturb positions slightly
sys.positions .= map(r -> r + 0.25 * randv(), sys.positions)
maybe_rebuild(sys, master; method = :cells)
master_after = length(master.pairs)
(master_before, master_after)
```

## Exclusions and skins

- Every pair potential carries an `exclusions::Vector{Tuple{Int,Int}}`. The
  default implementation performs a membership lookup when building the
  per-potential neighbour list.
- The per-potential `skin` is added to each parameter's cutoff when deciding
  whether to accept a pair from the master list. Choose `skin` large enough to
  avoid frequent rebuilds but small enough to keep lists tight.

## Diagnostics

- Use `brute_force_pairs(sys, cutoff)` to validate neighbour construction.
- Inspect `master.max_disp2` to see how far particles travelled since the last
  rebuild.
- `master.nbuilds` counts how many times the list has been rebuilt—handy for
  profiling.

For a hands-on tour, revisit [Tutorial 2](../tutorials/pair_potentials.md) which
relies on the automated neighbour management inside `ForceField`.

# Spec: Module `Verlet.Neighbors`

Purpose: Construct Verlet neighbour lists, manage shared neighbour state for forcefields, and provide utilities for validating pair selection.

## Core types
- `struct NeighborPair{F<:AbstractPotentialPair, IntT<:Integer}` stores one `(i, j)` pair together with its per-type parameters `pair::F`.
- `struct PotentialNeighborList{F, IntT}` wraps `neighbors::Vector{NeighborPair{F, IntT}}`; pair potentials create instances via `PotentialNeighborList(eltype(params.table))` so the stored element type matches their parameter table.
- `mutable struct ForceField{Layers}`
  - Fields: `layers::Layers`, `master::Any`, `master_method::Symbol`.
  - Constructor: `ForceField(layers; method=:cells)` initialises `master` to `nothing` (pair-potentials present) or `EmptyNeighborList()`.
  - Specialised `prepare_neighbors!(ff::ForceField, sys::System)` (extending `Verlet.Core.prepare_neighbors!`) ensures a `MasterNeighborList` exists, updates its `cutoff`/`skin`, rebuilds if needed, and then calls `build_neighbors_from_master!` for each pair potential in `layers`.
- `struct EmptyNeighborList end` acts as a sentinel when no pair potentials are present.

## Master neighbour list
- `mutable struct MasterNeighborList{D,T}`
  - Fields: `cutoff`, `cutoff2`, `skin`, `rskin2`, `pairs::Vector{SVector{2,Int}}`, `r0::Vector{SVector{D,T}}`, `max_disp2::T`, `nbuilds::Int`, `ncells::NTuple{D,Int}`, and reusable buffers `head::Vector{Int}`, `next::Vector{Int}` for the cell-linked build.
  - Constructors: `MasterNeighborList(sys::System; cutoff, skin)` and `MasterNeighborList(positions, box; cutoff, skin)`.
  - `rebuild!(nl, sys_or_positions; method=:cells, cutoff=nl.cutoff)` refreshes stored pairs using either the cell-linked sweep, a bruteforce fallback (`:bruteforce`), or an all-pairs builder (`:all_pairs`). Buffers are resized as needed while retaining allocations between calls.
  - `build_master_neighborlist!(nl, sys; r_verlet, method=:cells)` converts the requested Verlet radius to an internal cutoff (`max(0, r_verlet - nl.skin)` except for `:all_pairs`) before delegating to `rebuild!`.
- `rebuild_neighbors!(system::System, master::MasterNeighborList; kwargs...)` extends `Verlet.Core.rebuild_neighbors!` and simply calls `build_all_neighbors!(master, system.forcefield, system; kwargs...)`.
- `maybe_rebuild(system::System, master::MasterNeighborList; kwargs...)` (defined here via `Core.maybe_rebuild`) measures the maximum squared displacement of particles since `r0`; when it exceeds `(skin/2)^2`, the neighbour list is rebuilt.

## Building per-potential neighbours
- `build_neighbors_from_master!(pot::AbstractPairPotential, sys::System, master::MasterNeighborList)` clears `pot.neighborlist`, iterates `master.pairs`, applies minimum-image displacements, skips `is_excluded(pot, i, j)`, and keeps neighbours within `(rc + pot.skin)^2`.
- `is_excluded(pot::AbstractPairPotential, i, j)` performs a simple tuple-membership test on `pot.exclusions` (overridable for efficiency).
- `build_all_neighbors!(master::MasterNeighborList, ff::ForceField, sys::System; method=:cells)`
  - Collects pair potentials from `ff.layers`.
  - Ensures the master list covers the largest `(rc + skin)` amongst them.
  - Rebuilds each potential’s neighbour list from `master`.

## Supporting utilities
- `brute_force_pairs(sys::System, cutoff)` (and the `positions, box, cutoff` variant) returns all `(i, j)` pairs within the squared cutoff using the minimum-image convention; useful for validation.
- `_build_pairs_cells!`, `_build_pairs_bruteforce!`, and `_build_pairs_allpairs!` (internal) implement the various pair-generation strategies.
- Distance helpers `displacement` and `distance2_minimum_image` provide allocation-free minimum-image calculations for `SVector` inputs.

## Cell-grid helper
- `struct CellGrid{D,IT,T}`
  - Fields: `L`, `cell_size`, `dims`, `heads`, `next`.
  - Built via `build_cellgrid(R, box; cell_size)` which chooses a uniform grid no finer than the requested `cell_size` and immediately calls `rebin!`.
  - `rebin!(grid, R, box)` resets `heads`/`next`, bins positions into periodic cells, and can be reused between builds.

## Behaviour & invariants
- `master.pairs` always stores unique `(i < j)` pairs.
- `master.r0` matches the particle positions when the list was last rebuilt; `max_disp2` records the largest squared displacement seen since then.
- Per-potential neighbour lists are cleared before reuse so capacity is retained but contents are refreshed every rebuild.
- The `ForceField` keeps the most recent `MasterNeighborList` instance in `ff.master`; users may supply a pre-built master list when calling `build_all_neighbors!` directly.

## Example
```julia
using Verlet, StaticArrays

box = CubicBox(8.0)
R = [SVector{3}(randn(), randn(), randn()) for _ in 1:16]; wrap_positions!(R, box)
sys = System(R, fill(SVector{3}(0.0, 0.0, 0.0), 16), fill(SVector{3}(0.0, 0.0, 0.0), 16), ones(16), box,
             ones(Int, 16), Dict(1 => :A))

params = PairTable(fill(LJPair(1.0, 1.0, 2.5), (1, 1)))
lj = Verlet.Potentials.LennardJones(params, Tuple{Int,Int}[], 0.3)
ff = ForceField((lj,))
master = MasterNeighborList(sys; cutoff=2.5, skin=0.3)
build_all_neighbors!(master, ff, sys)
```

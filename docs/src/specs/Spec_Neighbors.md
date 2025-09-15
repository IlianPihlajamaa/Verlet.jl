# Spec: Module `Verlet.Neighbors`

Purpose: Efficient neighbor list construction and orchestration for force evaluation.

## Types

- `struct NeighborPair{F<:AbstractPotentialPair, IntT<:Integer}`
  - `i::IntT`, `j::IntT`, `pair::F` (per-pair parameters).
- `const PotentialNeighborList{F} = StructArray{NeighborPair{F,T_Int}} where {F<:AbstractPotentialPair}`
  - Per-potential neighbor container (StructArray-backed for SoA layout).
- `struct MasterNeighborEntry{T_Float,T_Int<:Integer}`
  - `i::T_Int`, `j::T_Int`, `r2::T_Float` (squared distance with minimum-image).
- `mutable struct MasterNeighborList{T<:Number} <: AbstractNeighborList`
  - `skin::T`, `entries::Vector{MasterNeighborEntry}`.
  - Constructor: `MasterNeighborList(skin; sizehint=1000)`.

## Building neighbors

- `build_master_neighborlist!(master, positions, box; r_verlet, method=:cells)`
  - Purpose: Fill `master.entries` with unique pairs `(i<j)` within cutoff `r_verlet` (plus algorithmic skin as needed).
  - Methods:
    - `:cells` (default): cell-list algorithm, expected `O(N)`.
    - `:bruteforce`: all `i<j` pairs, `O(N^2)`; useful for small N or debugging.
    - `:all_pairs`: include all pairs regardless of cutoff (testing/no-cutoff potentials).
  - Semantics: `r2` is computed with minimum-image in the provided `box`.
- `build_cellgrid(...)`, `rebin!(...)`
  - Utilities for the cell-list implementation (internal surface area; caller rarely uses directly).

## ForceField integration

- `build_neighbors_from_master!(pot::AbstractPairPotential, sys::System, master::MasterNeighborList)`
  - Populates `pot.neighborlist` using `master.entries` and per-type pair parameters.
  - Includes a pair `(i,j)` if not excluded and `entry.r2 < (p.rc + pot.skin)^2`.
- `build_all_neighbors!(master, ff::ForceField, sys::System; method=:cells)`
  - Computes a master list with `r_verlet = maximum(max(p.rc) + pot.skin for pot in ff.layers)` (over all layer parameter tables), then builds per-layer lists.

## Exclusions

- `is_excluded(pot::AbstractPairPotential, i, j)`
  - Basic tuple membership check on `pot.exclusions` (callers can precompute or specialize for performance).

## Invariants

- `master.entries` contains unique `(i<j)` pairs; no self-pairs; `r2 ≥ 0`.
- Per-potential neighbor lists are empty!+push!-rebuilt; capacity may be retained.

## Performance

- Cell-list path is expected `O(N)` with good binning; bruteforce is `O(N^2)`.
- StructArray neighbors enable tight loops in potentials (good cache behavior).

## Example

```julia
using Verlet, StaticArrays
box = CubicBox(10.0)
R = [@SVector randn(3) for _ in 1:128]; wrap_positions!(R, box)
sys = System(R, fill(@SVector zeros(3), 128), fill(@SVector zeros(3), 128), ones(128), box, ones(Int,128), Dict(1=>:A))

lj = begin
  params = PairTable(fill(LJPair(1.0, 1.0, 2.5), (1,1)))
  Verlet.Potentials.LennardJones(params, Tuple{Int,Int}[], 0.5)
end
ff = ForceField((lj,))
master = MasterNeighborList(0.5)
build_all_neighbors!(master, ff, sys; method=:cells)
```

## Cell Grid Details

- `CellGrid{IT,T}` fields:
  - `L::T`: cubic box length; assumed consistent with positions’ units.
  - `cell_size::T`: effective uniform cell width actually used for binning.
  - `dims::NTuple{3,IT}`: `(nx,ny,nz)`; each ≥ 1.
  - `heads::Vector{IT}`: length `nx*ny*nz`; head index per cell (`0` sentinel = empty).
  - `next::Vector{IT}`: length `N`; intrusive linked-list “next” pointer per particle (`0` = end).
- Build and rebin:
  - `build_cellgrid(R, box; cell_size)` constructs a grid sized so that the effective width `cs_eff = L / floor(L/cell_size) ≥ cell_size` and immediately `rebin!`s the positions.
  - `rebin!(grid, R, box)` resets `heads` and `next` and push-fronts each particle index `i` into its cell list using periodic mapping.
- Indexing and periodicity:
  - Positions are mapped to `[0, L)` via `x0 = x + 0.5L; x0 -= floor(x0 / L) * L`, then to cell indices `1..n` via `floor(x0/cell_size)+1` clamped to `[1,n]`.
  - Cell linear index: `((cz-1)*ny + (cy-1))*nx + cx`.
  - 27-neighbor stencil computed by `_neighbors_of_cell` using `mod1` periodic wrap per axis.

## Master List Algorithm

- Methods in `build_master_neighborlist!`:
  - `:cells` (default):
    - Ensure a compatible grid: if no grid provided, or `grid.cell_size < rlist` or `grid.L != box.L`, build a new grid with `cell_size=rlist`; otherwise `rebin!`.
    - For each cell: for each particle `i` in the cell’s list:
      - Same-cell pairs: traverse `j = next[i]` chain, test `r2 ≤ rlist2` using minimum-image, and `push!(i<j)`.
      - Neighbor cells: iterate the 27 neighbors; to avoid duplicates, only process neighbor cell `cc_idx` if `c_idx < cc_idx`, then traverse `heads[cc_idx]` list and test pairs.
    - After filling, `sort!` by `(i,j)` and `unique!` the entries.
  - `:bruteforce`: nested loops for `i<j`, compute `r2` with minimum-image and push if `r2 ≤ rlist2`.
  - `:all_pairs`: nested loops for `i<j`, always push with computed `r2` (no cutoff).
- Radii:
  - `rlist2 = (r_verlet + master_nl.skin)^2` for `:cells` and `:bruteforce`. `:all_pairs` ignores `rlist2`.

## Sorting and Uniqueness

- `:cells` path sorts and de-duplicates entries defensively, even though the neighbor sweep avoids double counting by `c_idx < cc_idx` and `i<j` ordering.
- `:bruteforce` and `:all_pairs` generate unique `(i<j)` pairs by construction and do not sort/unique.

## Grid Reuse Policy

- Pass an existing `grid` to `build_master_neighborlist!` to avoid reallocations.
- The grid is rebuilt if the requested list radius grows (`grid.cell_size < rlist`) or the box length changes; otherwise positions are re-binned in-place.

## Assumptions and Invariants

- Positions are 3D `SVector{3,T}`; the builder asserts `d=3`.
- `MasterNeighborEntry` stores `(i<j, r2≥0)` with `r2` computed under minimum-image in the provided `CubicBox`.
- `master.entries` is cleared (`resize!(…,0)`) before population.


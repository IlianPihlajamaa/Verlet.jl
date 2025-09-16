# Spec: Module `Verlet.Neighbors`

Purpose: Efficient neighbor list construction and orchestration for force evaluation.

## Types

- `struct NeighborPair{F<:AbstractPotentialPair, IntT<:Integer}`
  - `i::IntT`, `j::IntT`, `pair::F` (per-pair parameters).
- `const PotentialNeighborList{F} = StructArray{NeighborPair{F,T_Int}} where {F<:AbstractPotentialPair}`
  - Per-potential neighbor container (StructArray-backed for SoA layout).
- `mutable struct MasterNeighborList{D,T}`
  - Tracks `cutoff`, `skin`, and a `pairs::Vector{SVector{2,Int}}` holding unique `(i<j)` candidates within `cutoff + skin` using minimum-image distances.
  - Stores reference positions (`r0`), displacement bounds, the active cell grid shape `ncells`, and reusable `head`/`next` buffers for the cell build.
  - Constructors: `MasterNeighborList(sys; cutoff, skin)` or `MasterNeighborList(positions, box; cutoff, skin)`.
- Helper: `brute_force_pairs(sys_or_positions, box_or_cutoff, cutoff)` returns the set of `(i,j)` pairs within the provided cutoff using an `O(N^2)` sweep (useful for validation).

## Building neighbors

- `build_master_neighborlist!(master, positions_or_sys, box_or_kwargs; r_verlet, method=:cells)`
  - Purpose: Refresh `master.pairs` with unique `(i<j)` candidates consistent with the requested Verlet radius `r_verlet`.
  - `master.cutoff` is updated to `max(0, r_verlet - master.skin)` and `master.r0` holds the reference positions used for future displacement checks.
  - Methods:
    - `:cells` (default): cell-list algorithm using the internal `ncells`, `head`, and `next` buffers, expected `O(N)`.
    - `:bruteforce`: `O(N^2)` sweep that applies the cutoff directly; helpful for validation.
    - `:all_pairs`: fills every `(i<j)` pair ignoring the cutoff (testing / debugging).
- `rebuild!(master, sys_or_positions, box; method=:cells, cutoff=master.cutoff)`
  - Lower-level entry point called by `build_master_neighborlist!`; accepts either a `System` or raw positions plus box.
- `build_cellgrid(...)`, `rebin!(...)`
- `build_cellgrid(...)`, `rebin!(...)`
  - Utilities for the cell-list implementation (internal surface area; caller rarely uses directly).

## ForceField integration

- `build_neighbors_from_master!(pot::AbstractPairPotential, sys::System, master::MasterNeighborList)`
  - Populates `pot.neighborlist` using `master.pairs` and per-type pair parameters.
  - Computes the current squared distance on the fly (minimum-image) and includes `(i,j)` when not excluded and `r2 < (p.rc + pot.skin)^2`.
- `build_all_neighbors!(master, ff::ForceField, sys::System; method=:cells)`
  - Computes a master list with `r_verlet = maximum(max(p.rc) + pot.skin for pot in ff.layers)` (over all layer parameter tables), then builds per-layer lists.

## Exclusions

- `is_excluded(pot::AbstractPairPotential, i, j)`
  - Basic tuple membership check on `pot.exclusions` (callers can precompute or specialize for performance).

## Invariants

- `master.pairs` contains unique `(i<j)` pairs; no self-pairs; the stored cutoff equals the requested Verlet radius minus `skin`.
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
master = MasterNeighborList(sys; cutoff=2.5, skin=0.5)
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

- `:cells` (default):
  - `_choose_cell_grid` selects a `(nx, ny, nz)` stencil so that each cell edge is at most `cutoff + skin`.
  - Particles are assigned to cells by fractional coordinates and stored in intrusive lists via the `head`/`next` buffers.
  - For each base cell, the algorithm visits the 27-neighbor stencil with wraparound and emits `(i, j)` with `i < j` whenever the minimum-image squared distance is ≤ `(cutoff + skin)^2`.
  - `ncells`, `head`, and `next` are kept on the list object and reused between builds; buffers resize only when the particle count or grid shape changes.
- `:bruteforce`: nested loops for `i<j`, computing minimum-image distances directly and respecting the same cutoff radius (expected `O(N^2)`).
- `:all_pairs`: nested loops for `i<j` that bypass the cutoff entirely; useful for debugging and consistency tests.

## Sorting and Uniqueness

- All methods generate unique `(i<j)` pairs by construction; no post-processing `sort!`/`unique!` is performed.

## Assumptions and Invariants

- Positions are 3D `SVector{3,T}`; the builder asserts `d=3`.
- Minimum-image distances under the provided `CubicBox` determine pair inclusion.
- `master.pairs` is cleared before each rebuild; buffers grow amortized.

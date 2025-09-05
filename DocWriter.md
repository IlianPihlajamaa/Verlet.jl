## 📚 Doc Writer Prompt
```

You are the **Doc Writer**. Your job is to generate user-facing documentation for new features outlined in the current Design iteration.

Your job:

* Add docstrings to new functions, structs, or modules (Julia triple-quoted style).
* Update or create a Markdown files in `docs/` with usage examples and documentation.
* It should produce a valid documentation site with Documenter.jl
* Ensure examples are runnable Julia code snippets, perhaps with @example blocks
* Highlight performance tips, units, and pitfalls.
* Your goal is to write a unified diff patch (with \\\diff fencing) showing exactly what changes to apply. Ensure the docs build.

Rules:

* If you need to read a file to insert docstrings or inspect the file tree, request it using bash commands: Reply only with the command and nothing else and you shall receive the result. E.g.:

```bash
cat src/filename.jl
```


* If you want to run any code reply only with a bash command:

```bash
julia -e 'import Pkg; Pkg.test()'
```
Will provide the output of the tests. Similarly you can test building the documentation, and inspect the code.

* Always assume your working directory is the root of the package.

# Design
Here’s a **drop-in replacement** for your `Design.md` with the requested edits folded in (half-list emphasis + benchmark guidance + rebuild/skin notes). You can paste this over the current file.

---

# DESIGN.md — Next Feature Plan (revised)

## Overview — **Cell-Linked Lists (Bins) for O(N) Neighbor Build** + **Half Neighbor Lists**

We have a working Verlet neighbor list and LJ+PBC kernels. Benchmarks show NL force evaluation overtakes brute force around **N ≈ 256** and scales to **\~46×** at **N = 8192** (ρ=1, rcut=2.5, skin=0.4). However, the **current O(N²) neighbor build** dominates runtime as N grows.
Next, we add a **cell-linked grid** to make neighbor-list construction **O(N)** and switch the production path to a **half list** (store each pair once), halving memory and removing a branch in the kernel.

Scope is tight and non-breaking:

* Add `CellGrid` for cubic periodic boxes.
* Provide `build_cellgrid`, `rebin!`, and `build_neighborlist_cells` (O(N) builder).
* Emit **half neighbor lists** (pair stored once with `j>i`).
* Keep existing APIs; the cell-based builder is an alternative to the current O(N²) builder.

---

## Public API

```julia
struct CellGrid{IT<:Integer, T<:Real}
    L::T                # box length (CubicBox-compatible)
    cell_size::T        # typically rlist = cutoff + skin
    dims::NTuple{3,IT}  # (nx, ny, nz)
    heads::Vector{IT}   # length = nx*ny*nz; head of linked list per cell (0 = empty)
    next::Vector{IT}    # length = N; next particle index in cell list (0 = end)
end

build_cellgrid(R::AbstractMatrix, box::CubicBox; cell_size::Real) -> CellGrid
rebin!(grid::CellGrid, R::AbstractMatrix, box::CubicBox) -> CellGrid

build_neighborlist_cells(R::AbstractMatrix, box::CubicBox;
                         cutoff::Real, skin::Real=0.3,
                         grid::Union{Nothing,CellGrid}=nothing) -> NeighborList
# Emits a **half list**: each pair appears once with j>i.
```

**Exports (add):**

```julia
export CellGrid, build_cellgrid, rebin!, build_neighborlist_cells
```

Existing exports for `NeighborList`, `build_neighborlist`, `maybe_rebuild!`, `max_displacement_since_build`, `wrap_positions!` remain unchanged.

---

## Data Structures

```julia
struct CellGrid{IT<:Integer, T<:Real}
    L::T
    cell_size::T                     # enforce cell_size ≤ rlist internally
    dims::NTuple{3,IT}               # ≥ (1,1,1)
    heads::Vector{IT}                # 1-based indices; 0 sentinel
    next::Vector{IT}
end
```

* Linked list per cell avoids allocations during binning.
* `dims = max(1, floor(Int, L/cell_size))` per axis; clamp so `cell_size ≤ rlist`.
* Neighbor search visits the **27 neighboring cells** with periodic wrap.

---

## Algorithms

### A) Build / Rebin Grid — O(N)

1. Compute `(nx,ny,nz)` from `L` and `cell_size`.
2. Zero `heads` (length `nx*ny*nz`) and `next` (length `N`).
3. For each particle `i`, map position to **\[0,L)** then to `(cx,cy,cz)`, linearize to cell id `c`, and push-front: `next[i]=heads[c]; heads[c]=i`.

### B) Build Half Neighbor List from Cells — O(N) at fixed density

* `rlist = cutoff + skin`, precompute `rlist²`.
* For each cell and particle `i` in it:

  * For each of the 27 neighbor cells:

    * For each particle `j` in that cell: **skip if `j ≤ i`** (enforce half list).
    * Compute minimum-image `Δ`; if `‖Δ‖² ≤ rlist²`, record `(i,j)`.
* Assemble CSR: first pass counts → `offsets` (N+1), second pass fills `pairs`.
* Store `ref_positions = copy(R)`.

### C) List-aware LJ Kernel (half list)

Update the NL LJ kernel to **assume each pair appears once (half list)**. Remove runtime `j>i` checks and accumulate both particles’ forces in a single visit:

```
for i in 1:N
  for idx in offsets[i]:(offsets[i+1]-1)
    j = pairs[idx]  # j > i guaranteed by builder
    Δ = R[i,:] - R[j,:]; minimum_image!(Δ, box)
    r2 = dot(Δ,Δ)
    if r2 ≤ cutoff^2
        invr2 = 1/r2
        s2 = (σ^2)*invr2
        s6 = s2^3
        fr_over_r = 24*ϵ*(2*s6^2 - s6)*invr2
        Fi += fr_over_r * Δ
        Fj -= fr_over_r * Δ
        U  += 4*ϵ*(s6^2 - s6) - (shift ? Uc : 0)
    end
  end
end
```

### D) Rebuild Policy

Keep `maybe_rebuild!` (half-skin rule) unchanged. With O(N) builds, users can safely reduce `skin` to improve fidelity without large rebuild penalties.

---

## Numerical Pitfalls

* **Cell size**: enforce `cell_size ≤ rlist` or you risk missing pairs across cells.
* **Tiny boxes / low dims**: `dims` can be `(1,1,1)` → degenerates to O(N²) inside one cell (correct but slower).
* **Precision**: use `Float64` for distance math and accumulators.
* **Geometry guard**: recommend `L > 2*(cutoff+skin)` to avoid ambiguous minimum-image shells.
* **Half-list invariants**: kernel must not double-count; builder guarantees `j>i`.

---

## Acceptance Tests

Add `test/test_cellgrid.jl` and include it from `test/runtests.jl`.

```julia
@testset "CellGrid build and rebin invariants" begin
    N, D = 100, 3
    L = 12.0
    box = CubicBox(L)
    R = rand(N, D) .* L .- (L/2)
    wrap_positions!(R, box)

    rcut, skin = 2.5, 0.4
    rlist = rcut + skin

    grid = build_cellgrid(R, box; cell_size=rlist)
    @test grid.L == L
    @test all(x -> x ≥ 1, grid.dims)
    @test length(grid.next) == N
    @test length(grid.heads) == prod(grid.dims)

    rebin!(grid, R, box)
    @test length(grid.next) == N
end

@testset "Cell-based neighbor list ≡ naive rlist filter (half list)" begin
    N, D = 128, 3
    L = 15.0
    box = CubicBox(L)
    R = rand(N, D) .* L .- (L/2)
    wrap_positions!(R, box)

    cutoff, skin = 2.2, 0.5
    nlist = build_neighborlist_cells(R, box; cutoff=cutoff, skin=skin)

    rlist2 = (cutoff + skin)^2
    ref_pairs = Vector{Tuple{Int,Int}}()
    for i in 1:N-1, j in i+1:N
        Δ = @view R[i,:] .- R[j,:]
        minimum_image!(Δ, box)
        if dot(Δ,Δ) ≤ rlist2 + 1e-12
            push!(ref_pairs, (i,j))
        end
    end
    sort!(ref_pairs)

    got = Vector{Tuple{Int,Int}}()
    for i in 1:N
        for idx in nlist.offsets[i]:(nlist.offsets[i+1]-1)
            j = nlist.pairs[idx]
            @test j > i
            push!(got, (i,j))
        end
    end
    sort!(got)
    @test got == ref_pairs
end

@testset "LJ with half list from cells matches brute force" begin
    N, D = 64, 3
    L = 12.0
    box = CubicBox(L)
    R = rand(N, D) .* L .- (L/2)
    wrap_positions!(R, box)

    ϵ, σ, cutoff, skin = 1.0, 1.0, 2.5, 0.4
    nlist = build_neighborlist_cells(R, box; cutoff=cutoff, skin=skin)

    F_ref, U_ref = lj_forces(R, box; ϵ=ϵ, σ=σ, rcut=cutoff, return_potential=true)
    F_nl,  U_nl  = lj_forces(R, box, nlist; ϵ=ϵ, σ=σ, shift=false, return_potential=true)

    @test isapprox(F_nl, F_ref; atol=1e-10, rtol=1e-10)
    @test isapprox(U_nl,  U_ref; atol=1e-10, rtol=1e-10)
end

@testset "Rebin + maybe_rebuild! integration" begin
    N, D = 80, 3
    L = 14.0
    box = CubicBox(L)
    R = rand(N, D) .* L .- (L/2)
    wrap_positions!(R, box)

    cutoff, skin = 2.5, 0.4
    grid = build_cellgrid(R, box; cell_size=cutoff+skin)
    nlist = build_neighborlist_cells(R, box; cutoff=cutoff, skin=skin, grid=grid)

    R .+= 0.1 .* randn(N, D); wrap_positions!(R, box)
    @test maybe_rebuild!(nlist, R, box) == false

    R .+= 0.35 .* randn(N, D); wrap_positions!(R, box)
    if maybe_rebuild!(nlist, R, box)
        rebin!(grid, R, box)
        nlist2 = build_neighborlist_cells(R, box; cutoff=cutoff, skin=skin, grid=grid)
        @test length(nlist2.pairs) > 0
    end
end
```

---

## Benchmark Guidance (documentation, not a test)

* At **ρ≈1, D=3, rcut=2.5, skin=0.4**:

  * NL force evaluation overtakes O(N²) at **N≈256**.
  * Speedup grows with N: \~3× (512), \~6× (1024), \~23× (4096), \~46× (8192).
* With **O(N) builds**, you may reduce **skin** (e.g., 0.2–0.3) to improve force accuracy without paying a quadratic rebuild cost.

---

## Implementation Notes

* Files:

  * `src/cellgrid.jl`: `CellGrid`, `build_cellgrid`, `rebin!`, helpers (`cell_index`, periodic wrap).
  * `src/neighborlist_cells.jl`: `build_neighborlist_cells` (emits **half list**).
  * `include` both from `src/Verlet.jl`; add exports listed above.
* Ensure **half list** throughout: the cell builder guarantees `j>i`; update the NL LJ kernel to **not** branch on `j>i`.
* Keep distance math and accumulators in `Float64`.
* Enforce `cell_size ≤ rlist` inside `build_cellgrid` to preserve correctness.

---

## Task for Implementer (small & focused)

1. Implement `CellGrid`, `build_cellgrid`, and `rebin!` (O(N) binning) in `src/cellgrid.jl`.
2. Implement `build_neighborlist_cells` that:

   * Uses 27-cell sweep with periodic wrap,
   * Emits a **half list** with `j>i`,
   * Returns a valid CSR `NeighborList` with `ref_positions = copy(R)`.
3. Update `src/Verlet.jl` to include new files and exports.
4. Add `test/test_cellgrid.jl` with the Acceptance Tests and include from `test/runtests.jl`.
5. (Optional) Add `bench/bench_cellgrid_build.jl` to compare O(N²) vs O(N) builds at N∈{256,512,1024,2048}.

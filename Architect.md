## 🧱 Architect Prompt

````
You are the **MD Architect** working in the root directory of a Julia package Verlet (a molecular dynamics framework).
The end goal is to provide a fully featured MD engine in julia. It should be modern, extensible, flexible, performant and easy to use.
It must have modern features to differentiate itself from lammps/gromacs/etc. You can use other Julia packages if required/convenient. 
Also ensure that the code remains well structured.

Your job is to inspect the current state of the project and:
- Plan the next feature or fix.
- Edit the DESIGN document in Markdown with:
  - Overview (what & why)
  - Public API (functions, structs, modules, exports)
  - Data Structures (fields, types, mutability)
  - Algorithms (step-by-step, with pseudocode if needed)
  - Numerical Pitfalls (precision, stability issues)
  - Acceptance Tests (exact Julia `@test` conditions to add to `test/runtests.jl`)

Rules:
- Do not write actual implementation code.
- If you need to inspect existing files before designing, output only a shell command such as:
```bash
cat src/integrators.jl
````

* If you want to run code, you can request the output of shell commands by replying just the command and nothing else. e.g. reply by

```bash
julia -e 'import Pkg; Pkg.test()'
```
to see if the tests run.  

* Always assume your working directory is the root of the package.

* You can also inspect the current state of the package by visiting https://github.com/IlianPihlajamaa/Verlet.jl

* Once you have designed a good next step. Write a short task description in the DESIGN for the implementer, who will implement your vision. Ensure that it
is a small and resonable task. 

# Current DESIGN.md
# DESIGN.md — Next Feature Plan

## Overview — Add **Verlet Neighbor List** (O(N) average) for LJ + PBC

With a naive O(N²) Lennard–Jones force kernel and periodic cubic box now in scope, the next high-impact, contained step is a **Verlet neighbor list** (a.k.a. Verlet list) with a **skin** distance. This reduces the per-step force cost to \~O(N) for fixed density, enabling hundreds to thousands of particles at interactive speeds while preserving identical physics.

Scope is deliberately tight and orthogonal:

* A **`NeighborList`** type storing neighbors per particle for a given **cutoff** + **skin**.
* Builders: `build_neighborlist` (full rebuild) and `maybe_rebuild!` (displacement-based heuristic).
* A **list-aware LJ kernel**: `lj_forces(R, box, nlist; ...)` that iterates only over neighbors.
* Light **position wrapping** utility to keep diagnostics sane.
* Metadata and invariants for correctness & reproducibility.

No spatial hashing/cell lists yet; that can be a subsequent optimization once the API lands.

---

## Public API

```julia
# Types
struct NeighborList{IT<:Integer, T<:Real}
    cutoff::T       # physical cutoff used for pairs (rcut)
    skin::T         # additional shell beyond cutoff for list validity
    pairs::Vector{IT}      # flattened neighbor indices (CSR-style)
    offsets::Vector{IT}    # length N+1; neighbors for i are pairs[offsets[i]:(offsets[i+1]-1)]
    ref_positions::Matrix{T}  # copy of positions at last (re)build (N×D)
end
```

```julia
# Construction & maintenance
build_neighborlist(R::AbstractMatrix, box::CubicBox;
                   cutoff::Real, skin::Real=0.3) -> NeighborList

maybe_rebuild!(nlist::NeighborList, R::AbstractMatrix, box::CubicBox) -> Bool
# Returns true if a rebuild was performed (max displacement ≥ skin/2), false otherwise.

max_displacement_since_build(nlist::NeighborList, R::AbstractMatrix, box::CubicBox) -> Float64
```

```julia
# Force kernels
lj_forces(R::AbstractMatrix, box::CubicBox, nlist::NeighborList;
          ϵ::Real=1.0, σ::Real=1.0, shift::Bool=false, return_potential::Bool=false)

# Overload without nlist uses O(N²) reference implementation (already planned/added).
```

```julia
# Utilities
wrap_positions!(R::AbstractMatrix, box::CubicBox) -> R  # Wrap each coordinate into (-L/2, L/2]
```

**Exports (suggested):**

```julia
export NeighborList, build_neighborlist, maybe_rebuild!, max_displacement_since_build,
       wrap_positions!
```

(Keep `lj_forces` exported as already planned; method with `NeighborList` is additional multiple dispatch.)

---

## Data Structures

```julia
struct NeighborList{IT<:Integer, T<:Real}
    cutoff::T
    skin::T
    pairs::Vector{IT}      # concatenated neighbor indices (1-based)
    offsets::Vector{IT}    # CSR offsets; length = N+1
    ref_positions::Matrix{T}  # N×D snapshot at last build (Float64 default)
end
```

* **CSR layout** provides cache-friendly iteration and compact storage.
* `ref_positions` is used to track the **maximum particle displacement** since build (via minimum-image deltas), enabling cheap rebuild decisions.
* Use `Float64` internally for distances even if `R` is `Float32`, to keep numerical stability.

---

## Algorithms

### A. Building the neighbor list (naive O(N²) pass, still fast and simple)

Inputs: `R (N×D)`, `box`, `cutoff`, `skin`. Let `rlist = cutoff + skin`.

1. Precompute `rlist²`.
2. Allocate temporary `Vector{Vector{Int}} neighs(N)`.
3. For `i=1..N-1`, `j=i+1..N`:

   * `Δ = R[i,:] - R[j,:]`; `minimum_image!(Δ, box)`.
   * If `dot(Δ,Δ) ≤ rlist²`: push `j` into `neighs[i]` and `i` into `neighs[j]`.
4. Convert `neighs` to CSR:

   * `offsets[1]=1`, `offsets[i+1]=offsets[i]+length(neighs[i])`.
   * Write neighbors into `pairs`.
5. Store `ref_positions = copy(R)`.

This is a **one-time** or **infrequent** O(N²) step; per-step forces will be O(total neighbors) ≈ O(N).

### B. Rebuild heuristic (classic half-skin rule)

* Track **maximum** minimum-image displacement per particle since `ref_positions`.
* If `max_disp ≥ skin/2`, **rebuild** the list (`build_neighborlist`) and update `ref_positions`.
* Else keep using the list.

`maybe_rebuild!` implements this.

### C. List-aware LJ kernel

For each particle `i`:

```
Fi = 0
for j in neighbors(i):  # CSR via offsets/pairs
    Δ = R[i,:] - R[j,:]; minimum_image!(Δ, box)
    r2 = dot(Δ,Δ)
    if r2 <= cutoff^2
        invr2 = 1/r2
        s2 = (σ^2) * invr2
        s6 = s2^3
        fr_over_r = 24*ϵ*(2*s6^2 - s6)*invr2
        Fi += fr_over_r * Δ
        Fj -= fr_over_r * Δ  # enforce Newton-3rd; do once (i<j policy)
        U += 4*ϵ*(s6^2 - s6) - (shift ? Uc : 0)
    end
```

* Use **i\<j** policy when walking neighbors to avoid double-counting; for CSR we can either:

  * Build **symmetric** neighbor lists but apply interaction only when `i<j`, or
  * Build **half lists** (store each pair once). For simplicity & predictable cache patterns, we’ll use **symmetric lists + `if j>i` guard** in the first version.

### D. Position wrapping

`wrap_positions!` maps each coordinate of each particle to `(-L/2, L/2]`. This is **not required** by the kernel (since we minimum-image displacements), but improves debuggability and ensures `max_displacement_since_build` remains well-behaved.

---

## Numerical Pitfalls

* **Skin too small → frequent rebuilds** and potential missed pairs if particles cross the shell before rebuild. Default `skin=0.3σ` is a reasonable start; users should tune by temperature/density.
* **Precision**: Accumulate forces and energy in `Float64`. Displacement tracking must use minimum-image deltas to avoid spuriously large displacements across boundaries.
* **List symmetry & Newton’s third law**: If symmetric lists are used, ensure pair is applied once (e.g., `j>i`) or apply twice with half force—prefer the former.
* **Cutoff consistency**: The list uses `rlist = cutoff + skin`, but the **force** uses `cutoff` exactly. If `shift=true`, compute `Uc = U(cutoff)` (not at `rlist`).
* **Degenerate boxes**: Require `L > 2*(cutoff+skin)`; otherwise minimum-image is ambiguous and the neighbor shell spans the box.

---

## Acceptance Tests

Add these to `test/test_neighborlist.jl` and include from `runtests.jl`.

```julia
using Test, LinearAlgebra, Verlet

@testset "NeighborList construction and invariants" begin
    N, D = 10, 3
    L = 8.0
    box = CubicBox(L)
    R = randn(N, D)
    wrap_positions!(R, box)

    cutoff = 2.5
    skin   = 0.4
    nlist = build_neighborlist(R, box; cutoff=cutoff, skin=skin)

    # Offsets form valid CSR
    @test length(nlist.offsets) == N + 1
    @test first(nlist.offsets) == 1
    @test last(nlist.offsets) == length(nlist.pairs) + 1

    # All neighbor pairs within rlist
    rlist2 = (cutoff + skin)^2
    for i in 1:N
        for idx in nlist.offsets[i]:(nlist.offsets[i+1]-1)
            j = nlist.pairs[idx]
            Δ = @view R[i,:] .- R[j,:]
            minimum_image!(Δ, box)
            @test dot(Δ,Δ) <= rlist2 + 1e-12
        end
    end
end

@testset "maybe_rebuild! triggers on half-skin" begin
    N, D = 5, 3
    L = 10.0
    box = CubicBox(L)
    R = zeros(N, D)
    cutoff, skin = 2.0, 0.6
    nlist = build_neighborlist(R, box; cutoff=cutoff, skin=skin)

    # Move one particle by < skin/2 : no rebuild
    R[1,1] += 0.25  # < 0.3
    @test maybe_rebuild!(nlist, R, box) == false

    # Cross skin/2 : triggers rebuild
    R[1,1] += 0.1  # now 0.35 > 0.3
    @test maybe_rebuild!(nlist, R, box) == true
end

@testset "LJ forces with NeighborList match O(N²) reference" begin
    N, D = 24, 3
    L = 12.0
    box = CubicBox(L)
    R = rand(N, D) .* (L/2) .- (L/4)
    wrap_positions!(R, box)

    ϵ, σ, cutoff, skin = 1.0, 1.0, 2.5, 0.4
    # Build list with rlist = cutoff + skin
    nlist = build_neighborlist(R, box; cutoff=cutoff, skin=skin)

    # Reference (brute-force) and list-aware comparisons
    F_ref, U_ref = lj_forces(R, box; ϵ=ϵ, σ=σ, rcut=cutoff, return_potential=true)
    F_nl,  U_nl  = lj_forces(R, box, nlist; ϵ=ϵ, σ=σ, shift=false, return_potential=true)

    @test isapprox(F_nl, F_ref; atol=1e-10, rtol=1e-10)
    @test isapprox(U_nl,  U_ref; atol=1e-10, rtol=1e-10)
end

@testset "NeighborList stays valid under small motion; rebuild changes content" begin
    N, D = 16, 3
    L = 10.0
    box = CubicBox(L)
    R = rand(N, D) .* L .- (L/2)
    cutoff, skin = 2.2, 0.6
    nlist = build_neighborlist(R, box; cutoff=cutoff, skin=skin)

    pairs_before = copy(nlist.pairs)
    offsets_before = copy(nlist.offsets)

    # Small random motion below half-skin → no rebuild
    R .+= 0.1 .* randn(N, D)
    wrap_positions!(R, box)
    @test maybe_rebuild!(nlist, R, box) == false
    @test nlist.pairs == pairs_before
    @test nlist.offsets == offsets_before

    # Larger motion → rebuild likely alters neighbor graph
    R .+= 0.4 .* randn(N, D) # push beyond skin/2
    wrap_positions!(R, box)
    did = maybe_rebuild!(nlist, R, box)
    @test did == true
end

@testset "wrap_positions! maps coords to (-L/2, L/2]" begin
    box = CubicBox(5.0)
    R = [ 3.1 -2.6; -4.9 4.9 ]
    wrap_positions!(R, box)
    @test all(-box.L/2 .< R) && all(R .<= box.L/2)
end

@testset "Energy shift at cutoff is preserved with NeighborList" begin
    L = 30.0; box = CubicBox(L)
    R = [0.0 0.0 0.0; 2.0 0.0 0.0]
    cutoff, skin = 2.5, 0.4
    nlist = build_neighborlist(R, box; cutoff=cutoff, skin=skin)
    F1, U1 = lj_forces(R, box, nlist; rcut=cutoff, shift=true, return_potential=true)
    R2 = [0.0 0.0 0.0; 2.51 0.0 0.0]
    nlist2 = build_neighborlist(R2, box; cutoff=cutoff, skin=skin)
    F2, U2 = lj_forces(R2, box, nlist2; rcut=cutoff, shift=true, return_potential=true)
    @test isapprox(norm(F2), 0.0; atol=1e-12)
    @test U2 ≈ 0.0 atol=1e-10
    @test U1 < 0.0
end
```

---

## Implementation Notes (for the implementer)

* Put neighbor-list code in `src/neighborlist.jl`; `include` from `src/Verlet.jl`.
* Default `skin=0.3` (in units of `σ`, but **do not** couple to `σ` in code; just a numeric default).
* Keep allocations minimal in the LJ kernel by reusing a `Δ::Vector{Float64}` scratch.
* In `maybe_rebuild!`, compute per-particle displacement as `norm(minimum_image!(R[i,:] - ref_positions[i,:], box))` and track the maximum.
* Add `wrap_positions!` to `src/box.jl` (or alongside `CubicBox`).

---

## Task for Implementer (small & focused)

1. **Add types & functions**:

   * `NeighborList`, `build_neighborlist`, `maybe_rebuild!`, `max_displacement_since_build`, `wrap_positions!`.
   * Extend `lj_forces` with a method that accepts `(R, box, nlist; ...)` and iterates neighbors (respecting `j>i` to avoid double counting).
2. **Exports**: as suggested above.
3. **Tests**: create `test/test_neighborlist.jl` with the acceptance tests and `include("test_neighborlist.jl")` from `test/runtests.jl`.
4. **Docs**: a short section “Speeding up LJ with a Neighbor List” showing:

   ```julia
   box = CubicBox(20.0)
   nlist = build_neighborlist(R, box; cutoff=2.5, skin=0.4)
   maybe_rebuild!(nlist, R, box) && @info "rebuilt"
   F = lj_forces(R, box, nlist; rcut=2.5)
   ```
5. **(Optional)**: add a micro-benchmark in `bench/` comparing O(N²) vs neighbor list at N≈128–512.

Once this lands and `lj_forces` (brute-force) is in place, we can consider **cell lists** for O(N) builds and **time-correlated rebuild heuristics** as the next iteration.

# Implemented!

benchmarking:

#!/usr/bin/env julia
# Benchmark neighbor-list vs brute-force LJ with compilation & rebuild amortized
#
# Run with:
#   julia --project -O3 bench/bench_neighborlist.jl
#
# Notes:
# - Performs a warm-up phase to compile kernels before timing.
# - Benchmarks per-step cost over a trajectory with small random motions.
# - Uses identical motion for both methods to ensure fair comparison.

using Random, Statistics, Printf
using Verlet

struct BenchResult
    N::Int
    steps::Int
    warmup::Int
    brute_mean::Float64
    nlist_mean::Float64
    speedup::Float64
    rebuilds::Int
end

# Generate a reproducible set of displacements for all steps (warmup + measured).
function make_displacements(N, D, total_steps; disp=0.02, rng=Random.default_rng())
    disps = [disp .* randn(rng, N, D) for _ in 1:total_steps]
    return disps
end

function run_bruteforce(R0, box; disps::Vector, rcut::Real)
    R = copy(R0)
    # Warm compile
    Verlet.lj_forces(R, box; rcut=rcut)
    GC.gc()
    t = @elapsed begin
        for Δ in disps
            @inbounds R .+= Δ
            Verlet.wrap_positions!(R, box)
            Verlet.lj_forces(R, box; rcut=rcut)
        end
    end
    return t / length(disps)
end

function run_nlist(R0, box; disps::Vector, rcut::Real, skin::Real)
    R = copy(R0)
    nlist = Verlet.build_neighborlist(R, box; cutoff=rcut, skin=skin)
    # Warm compile
    Verlet.maybe_rebuild!(nlist, R, box)
    Verlet.lj_forces(R, box, nlist; rcut=rcut)
    GC.gc()
    rebuilds = 0
    t = @elapsed begin
        for Δ in disps
            @inbounds R .+= Δ
            Verlet.wrap_positions!(R, box)
            rebuilds += Verlet.maybe_rebuild!(nlist, R, box)
            Verlet.lj_forces(R, box, nlist; rcut=rcut)
        end
    end
    return (t / length(disps), rebuilds)
end

function bench_once(; N=256, D=3, L=20.0, rcut=2.5, skin=0.4,
                     disp=0.02, warmup=100, steps=400, seed=42)
    rng = MersenneTwister(seed)
    box = Verlet.CubicBox(L)
    R0 = randn(rng, N, D)
    Verlet.wrap_positions!(R0, box)
    total = warmup + steps
    disps = make_displacements(N, D, total; disp=disp, rng=rng)
    # Split into warmup + measured segments
    disps_warm = disps[1:warmup]
    disps_meas = disps[warmup+1:end]

    # Run separate warmups to compile both paths
    _ = run_bruteforce(R0, box; disps=disps_warm, rcut=rcut)
    _ = run_nlist(R0, box; disps=disps_warm, rcut=rcut, skin=skin)

    brute_mean = run_bruteforce(R0, box; disps=disps_meas, rcut=rcut)
    nlist_mean, rebuilds = run_nlist(R0, box; disps=disps_meas, rcut=rcut, skin=skin)

    return BenchResult(N, steps, warmup, brute_mean, nlist_mean, brute_mean / nlist_mean, rebuilds)
end

function bench_sweep(Ns = (64, 128, 256, 512, 1024);
                     steps=800, warmup=200, L=20.0, rcut=2.5, skin=0.4, disp=0.02, seed=42)
    results = BenchResult[]
    for N in Ns
        res = bench_once(N=N, L=L, rcut=rcut, skin=skin, disp=disp, steps=steps, warmup=warmup, seed=seed)
        push!(results, res)
    end
    return results
end

function print_results(results::Vector{BenchResult})
    @printf("%-6s %-8s %-8s %-12s %-12s %-8s %-9s\n",
            "N", "steps", "warmup", "brute/step", "nlist/step", "speedup", "rebuilds")
    for r in results
        @printf("%-6d %-8d %-8d %-12.6e %-12.6e %-8.2f %-9d\n",
                r.N, r.steps, r.warmup, r.brute_mean, r.nlist_mean, r.speedup, r.rebuilds)
    end
end

# Defaults are reasonable; tweak via environment variables if desired.
Ns = (64, 128, 256, 512, 1024)
results = bench_sweep(Ns;
                        steps=parse(Int, get(ENV, "STEPS", "800")),
                        warmup=parse(Int, get(ENV, "WARMUP", "200")),
                        L=parse(Float64, get(ENV, "L", "20.0")),
                        rcut=parse(Float64, get(ENV, "RCUT", "2.5")),
                        skin=parse(Float64, get(ENV, "SKIN", "0.4")),
                        disp=parse(Float64, get(ENV, "DISP", "0.02")),
                        seed=parse(Int, get(ENV, "SEED", "42")))
print_results(results)
N      steps    warmup   brute/step   nlist/step   speedup  rebuilds 
64     800      200      4.987998e-05 4.700798e-05 1.06     86       
128    800      200      2.017386e-04 1.767412e-04 1.14     97       
256    800      200      7.930778e-04 7.182588e-04 1.10     108      
512    800      200      3.161800e-03 3.009996e-03 1.05     118      
1024   800      200      1.261396e-02 1.250826e-02 1.01     127  
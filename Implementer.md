## 💻 Implementer Prompt
```

You are the **Implementer**. Your task is to implement the DESIGN provided by the Architect in the Julia package Verlet.

Your job:

* Modify relevant files under `src/`.
* Add or update tests in `test/runtests.jl`.
* Ensure the code type-checks and runs under Julia.
* Return a unified diff patch (with \`\`\`diff fencing) showing exactly what changes to apply. 

Rules:

* Do not modify unrelated code.
* Do not relax tests; make the code satisfy the acceptance tests.
* If you need to inspect a file, request it using only a shell code block with `cat`.
* If you want to run tests, request it with:

```bash
julia -e 'import Pkg; Pkg.test()'
```

* Always assume your working directory is the root of the package.
* If you want to run code, you can request the output of shell commands by replying just the command and nothing else. e.g. reply by

```bash
julia -e 'import Pkg; Pkg.test()'
```
to see if the tests run. You can also inspect the current state of the code by visiting https://github.com/IlianPihlajamaa/Verlet.jl

# Current DESIGN.md

# DESIGN.md — Next Feature Plan

## Overview — Add Periodic Box + Lennard–Jones Forces (O(N²))

To run physically meaningful MD, we need **inter-particle forces** and **boundary conditions**.
Next step: implement a **cubic periodic box** with **minimum-image** convention and a standard **Lennard–Jones (LJ) pair potential** with optional cutoff. This unlocks simple fluids, sanity checks for energy conservation, and a base for future optimizations (neighbor lists, cell lists, PME, etc.).

Scope is intentionally minimal and self-contained:

* **Data type** for cubic box.
* **Utility** to apply minimum-image displacement.
* **Force provider** `lj_forces` compatible with the existing `forces(r; return_potential=true)` convention used by `potential_energy`.

No neighbor list yet (keep it small and clear).

---

## Public API

### New types

* `struct CubicBox{T<:Real}`

  * `L::T` — box length (same in all dimensions).

### New functions

* `minimum_image!(Δ::AbstractVector, box::CubicBox)`

  * In-place minimum-image wrap of a displacement vector `Δ` so each component lies in `(-L/2, L/2]`.

* `lj_forces(positions::AbstractMatrix, box::CubicBox;
             ϵ::Real=1.0, σ::Real=1.0, rcut::Real=Inf,
             shift::Bool=false, return_potential::Bool=false)`

  * Returns an `(N×D)` force matrix; if `return_potential=true`, returns `(F, U::Float64)`.
  * **Pair potential:**
    $U(r) = 4\varepsilon \left[(\sigma/r)^{12} - (\sigma/r)^6\right]$
  * **Pair force magnitude along $\hat r$:**
    $f(r) = 24\varepsilon \left[2(\sigma/r)^{12} - (\sigma/r)^6\right]/r$
  * **Cutoff:** include pairs with `r ≤ rcut`.
    If `shift=true`, add a constant so that $U(rcut)=0$ (energy continuity); **no force smoothing** (derivative discontinuity remains).

These are **not exported** unless you prefer; suggested exports:

```julia
export CubicBox, minimum_image!, lj_forces
```

---

## Data Structures

```julia
struct CubicBox{T<:Real}
    L::T  # box length
end
```

* Immutable, trivially small.
* Assumes simulation dimensionality `D` equals `size(positions, 2)`; we do not store `D` in the box.

**Displacement scratch:**

* Internally, `lj_forces` uses a small `Vector{Float64}` (length `D`) as a reusable scratch for `Δ`.

---

## Algorithms

### Minimum-image displacement

For each component $\Delta_k$ of displacement vector `Δ`:

```
half = box.L / 2
if Δ[k] >  half: Δ[k] -= box.L
if Δ[k] ≤ -half: Δ[k] += box.L
```

This maps to `(-L/2, L/2]`. Use `while` only if you ever expect |Δ|>L (shouldn’t occur for pairwise displacements formed via ri - rj).

### Lennard–Jones forces (naive O(N²))

Inputs: positions `R::(N×D)`, `box::CubicBox`, parameters `(ϵ, σ, rcut, shift)`.

Pseudocode:

```
function lj_forces(R, box; ϵ=1, σ=1, rcut=Inf, shift=false, return_potential=false)
    N, D = size(R)
    F = zeros(N, D)
    U = 0.0

    σ2 = σ^2
    rcut2 = rcut^2

    # Optional energy shift
    Uc = 0.0
    if shift && isfinite(rcut)
        s2c = σ2 / rcut2
        s6c = s2c^3
        Uc  = 4*ϵ*(s6c^2 - s6c)  # U(rcut)
    end

    Δ = zeros(D)  # scratch

    for i in 1:N-1
        ri = @view R[i, :]
        for j in i+1:N
            rj = @view R[j, :]
            Δ .= ri .- rj
            minimum_image!(Δ, box)
            r2 = dot(Δ, Δ)
            if r2 <= rcut2
                invr2 = 1 / r2
                s2 = σ2 * invr2
                s6 = s2^3
                # Force magnitude over r: fr_over_r = 24ϵ*(2*s6^2 - s6) * invr2
                fr_over_r = 24*ϵ*(2*s6^2 - s6) * invr2
                # Vector force on i: Fi += fr_over_r * Δ
                for k in 1:D
                    f = fr_over_r * Δ[k]
                    F[i,k] += f
                    F[j,k] -= f
                end
                if return_potential
                    U += 4*ϵ*(s6^2 - s6) - Uc
                end
            end
        end
    end
    return return_potential ? (F, U) : F
end
```

Notes:

* All math avoids `sqrt` by working with `r²`; only the **unit direction** enters through multiplication by `Δ` (equivalent to dividing by `r`).
* Newton’s third law is enforced by symmetric accumulation (`+` for `i`, `-` for `j`).

---

## Numerical Pitfalls

* **Singularity at r→0**: huge forces; users must avoid overlapping particles or add softening. We guard only via finite arithmetic.
* **Cutoff discontinuity**: With `shift=false`, potential and force both jump at `rcut`. With `shift=true`, potential is continuous but **force remains discontinuous**. This can introduce small energy drift—acceptable for this minimal feature.
* **Finite precision**: Use `Float64` for positions/velocities to reduce drift; consistent with current package.
* **Box length vs positions**: Users must keep positions wrapped or large drifts can produce |Δ|>L/2 that violate minimum-image assumptions for direct `ri - rj`. Here we explicitly minimum-image the displacement, so positions themselves need not be wrapped each step, but keeping them within box helps debug output.
* **Performance**: O(N²) scales poorly; neighbor lists and cell lists are the natural next step.

---

## Acceptance Tests

Add to `test/runtests.jl` (or a new `test/test_lj.jl` included by `runtests.jl`). Use `atol=1e-10` where needed.

```julia
using Test
using Verlet

@testset "CubicBox + minimum_image!" begin
    box = CubicBox(10.0)
    Δ = [ 6.0, -6.0,  0.1]  # in 3D
    minimum_image!(Δ, box)
    @test Δ[1] ==  6.0 - 10.0    # -> -4.0
    @test Δ[2] == -6.0 + 10.0    # ->  4.0
    @test isapprox(Δ[3], 0.1; atol=1e-12)
end

@testset "LJ two-particle (no PBC, analytic check)" begin
    # Place two particles distance r along x
    r  = 1.5
    σ  = 1.0
    ϵ  = 2.0
    R  = [0.0 0.0 0.0;
          r   0.0 0.0]
    box = CubicBox(100.0) # effectively no wrapping
    F, U = lj_forces(R, box; ϵ=ϵ, σ=σ, rcut=Inf, return_potential=true)

    # Analytic force magnitude
    s    = σ/r
    s6   = s^6
    fmag = 24*ϵ*(2*s6^2 - s6)/r  # along +x on particle 1
    @test isapprox(F[1,1],  fmag; atol=1e-10)
    @test isapprox(F[2,1], -fmag; atol=1e-10)
    @test isapprox(F[1,2], 0.0; atol=1e-12)
    @test isapprox(F[2,3], 0.0; atol=1e-12)

    # Analytic potential
    Uref = 4*ϵ*(s6^2 - s6)
    @test isapprox(U, Uref; atol=1e-10)
end

@testset "LJ minimum-image under PBC" begin
    # Put particles near opposite faces of a small box; distance should wrap
    L = 5.0
    box = CubicBox(L)
    R = [ -2.4   0.0  0.0;   # ~ -L/2 + 0.1
           2.4   0.0  0.0]   # ~  L/2 - 0.1
    F = lj_forces(R, box; ϵ=1.0, σ=1.0, rcut=Inf, return_potential=false)
    # Displacement should be 0.2 along x after wrapping; force pushes apart
    @test F[1,1] > 0.0
    @test F[2,1] < 0.0
    @test isapprox(F[1,1], -F[2,1]; atol=1e-12)
end

@testset "LJ cutoff and energy shift" begin
    box = CubicBox(50.0)
    R   = [0.0 0.0; 1.2 0.0]  # 2D for variety
    rcut = 1.25
    # With shift
    F1, U1 = lj_forces(R, box; rcut=rcut, shift=true, return_potential=true)
    # Move just beyond cutoff: forces zero, potential ~0 with shift=true
    R2  = [0.0 0.0; 1.26 0.0]
    F2, U2 = lj_forces(R2, box; rcut=rcut, shift=true, return_potential=true)
    @test isapprox(norm(F2), 0.0; atol=1e-12)
    @test U2 ≈ 0.0 atol=1e-10
    @test U1 < 0.0  # attractive well inside cutoff
end

@testset "Compatibility with velocity_verlet! API" begin
    # Use lj_forces as the `forces` function with potential return
    box = CubicBox(20.0)
    forces(r; return_potential=false) = return_potential ?
        lj_forces(r, box; return_potential=true) :
        lj_forces(r, box; return_potential=false)

    ps = ParticleSystem([0.9 0.0 0.0; 0.0 0.0 0.0],
                        zeros(2,3),
                        ones(2))
    E0 = kinetic_energy(ps) + potential_energy(ps, forces)
    velocity_verlet!(ps, forces, 1e-3)
    E1 = kinetic_energy(ps) + potential_energy(ps, forces)
    @test isfinite(E0) && isfinite(E1)
end
```

---

## Task for Implementer (small & focused)

1. **Add types & functions** to `src/Verlet.jl` (or split into `src/forces.jl` and `include` it from `Verlet.jl`):

   * `CubicBox{T}`, `minimum_image!`, `lj_forces` exactly per API above.
   * Keep `lj_forces` free-function (not tied to `ParticleSystem`) so it composes with current `potential_energy` convention.
2. **Export** `CubicBox, minimum_image!, lj_forces` (optional but recommended).
3. **Write tests** as specified (you can place them in `test/test_lj.jl` and include from `test/runtests.jl`).
4. **Docs**: add a short usage snippet to `docs/src/guide/forces.md` demonstrating `forces(r) = lj_forces(r, CubicBox(…))`.

That’s it. Once this lands, we can profile and then design a **Verlet neighbor list** with rebuild heuristics and vectorized kernels as the next incremental feature.

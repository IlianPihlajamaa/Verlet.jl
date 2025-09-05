## 💻 Implementer Prompt
```

You are the **Implementer**. Your task is to implement the DESIGN provided by the Architect in the Julia package Verlet.

Your job:

* Modify relevant files under `src/`.
* Add or update tests in `test/runtests.jl`.
* Ensure the code type-checks and runs under Julia.
* Return a unified diff patch (with \`\`\`diff fencing) showing exactly what changes to apply. 
* Check that all acceptance tests in the DESIGN are covered,ensure the tests run.

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
julia -e 'import Pkg; Pkg.activate("."); Pkg.test()'
```
to see if the tests run. You can also inspect the current state of the code by visiting https://github.com/IlianPihlajamaa/Verlet.jl

# Current DESIGN.md
Here’s a **drop-in replacement** for your `Design.md` with the requested edits folded in (half-list emphasis + benchmark guidance + rebuild/skin notes). You can paste this over the current file.

---
# DESIGN.md — Next Feature Plan

## Overview — **Holonomic Distance Constraints (SHAKE/RATTLE) + COM Drift Removal**

With neighbor lists, LJ + PBC, and an NVT (BAOAB) thermostat on deck, the next capability that unlocks **larger stable timesteps** and **molecular realism** is support for **holonomic distance constraints** (e.g., rigid bonds to H). We’ll implement classical **SHAKE** (position projection) and **RATTLE** (velocity projection) for **pairwise distance constraints**, plus a light utility to **remove center-of-mass (COM) drift**.

Benefits:

* Enable typical “rigid bonds to H” workflows and **increase `dt`** (e.g., from 1 fs → 2 fs equivalents in reduced units).
* Provide a clean path to rigid molecules (e.g., TIP3P water via distance-only constraints first; angle/SETTLE can come later).
* Integrate with existing integrators; provide a **constrained VV** driver and a **constraint-aware DoF** for temperature tools.

Scope (tight, non-breaking):

* New `constraints.jl` module with a minimal `DistanceConstraints` struct.
* `velocity_verlet_shake_rattle!` integrator (wrapper over current VV) that applies SHAKE/RATTLE.
* Update `degrees_of_freedom` to account for constraints (optionally COM removal).
* Utility: `remove_com_motion!` (mass-weighted).

---

## Public API

```julia
"""
    DistanceConstraints(pairs, lengths; tol=1e-8, maxiter=50, use_minimum_image=true)

Create a distance-constraint set:
- pairs::Vector{Tuple{Int,Int}}: constrained atom index pairs (1-based)
- lengths::Vector{Float64}: target distances (same units as positions)
- tol: max |constraint violation| tolerated for convergence (in length^2)
- maxiter: maximum SHAKE and RATTLE iterations per step
- use_minimum_image: if true, constraints use minimum-image displacement under PBC
"""
struct DistanceConstraints end  # concrete fields described below
DistanceConstraints(pairs, lengths; tol=1e-8, maxiter=50, use_minimum_image=true)

"""
    velocity_verlet_shake_rattle!(ps::ParticleSystem, forces, dt, cons::DistanceConstraints)

Advance one constrained step:
1) VV half-kick, drift
2) SHAKE position projection
3) Recompute forces
4) VV half-kick
5) RATTLE velocity projection
"""
velocity_verlet_shake_rattle!(ps, forces, dt, cons)

"""
    apply_shake!(ps::ParticleSystem, cons::DistanceConstraints, dt)

Project positions to satisfy all constraints (used inside integrators).
"""
apply_shake!(ps, cons, dt)

"""
    apply_rattle!(ps::ParticleSystem, cons::DistanceConstraints)

Project velocities to satisfy the differential constraints d/dt C_l = 0.
"""
apply_rattle!(ps, cons)

"""
    degrees_of_freedom(ps; constraints=nothing, remove_com=false) -> Int

Return the effective translational DoF given optional constraints and COM removal.
"""
degrees_of_freedom(ps; constraints=nothing, remove_com=false) -> Int

"""
    remove_com_motion!(ps; which=:velocity)

Remove center-of-mass motion:
- which = :velocity | :position | :both
"""
remove_com_motion!(ps; which=:velocity)
```

**Exports (add):**

```julia
export DistanceConstraints,
       velocity_verlet_shake_rattle!,
       apply_shake!, apply_rattle!,
       remove_com_motion!
# `degrees_of_freedom` is already exported; extend its method with kwargs.
```

**Files (new/updated):**

* `src/constraints.jl` — all constraint types & algorithms.
* `src/Verlet.jl` — `include("constraints.jl")`, export new symbols, extend `degrees_of_freedom`.

---

## Data Structures

```julia
struct DistanceConstraints
    i::Vector{Int}          # first atom index per constraint
    j::Vector{Int}          # second atom index per constraint
    r0::Vector{Float64}     # target distances per constraint
    tol::Float64            # convergence tolerance on |C_l|
    maxiter::Int            # maximum iterations
    use_minimum_image::Bool # apply minimum-image to displacement vectors
end
```

* **Mutability:** keep `DistanceConstraints` **immutable** (parameters). All corrections are applied to `ParticleSystem` in place.
* **Masses:** read from `ps.masses`.
* **Box / PBC:** assume an orthorhombic box already exists in the codebase; if not, use raw positions and set `use_minimum_image=false` (documented).

---

## Algorithms

### Constraints

Each distance constraint `l` ties atoms `i,j` with

```
C_l(r) = ||r_i - r_j||^2 - r0_l^2 = 0.
```

#### SHAKE (positions)

After the unconstrained drift, we correct positions with Lagrange multipliers `{λ_l}`:

For constraint `l` with `d_ij = r_i - r_j` (use minimum image if requested),

```
σ_l = (1/m_i + 1/m_j) * (d_ij ⋅ d_ij)
Δλ_l = - C_l(r) / (2 * σ_l)
Δr_i =  (Δλ_l / m_i) * d_ij
Δr_j = -(Δλ_l / m_j) * d_ij
```

Apply Gauss–Seidel style updates **iteratively** over all constraints until:

```
max_l |C_l(r)| ≤ tol   or   iterations ≥ maxiter
```

**Notes:**

* Use current `d_ij` each sub-iteration (Gauss–Seidel converges better than Jacobi).
* If `σ_l` is extremely small (atoms nearly coincident), abort with a descriptive error.
* With PBC, compute `d_ij` via minimum image, but apply corrections to wrapped positions directly (positions should remain inside the primary cell; the constraint is topologically within a molecule, so this is physically consistent).

#### RATTLE (velocities)

Enforce the velocity constraint

```
d/dt C_l = 2 d_ij ⋅ (v_i - v_j) = 0.
```

Solve for `μ_l` and correct velocities:

```
τ_l = (1/m_i + 1/m_j) * (d_ij ⋅ d_ij)
μ_l = - (d_ij ⋅ (v_i - v_j)) / τ_l
Δv_i =  (μ_l / m_i) * d_ij
Δv_j = -(μ_l / m_j) * d_ij
```

Iterate over constraints until `max_l |d_ij ⋅ (v_i - v_j)| ≤ tol_v`, with `tol_v` tied to `tol` (e.g., `tol_v = tol^(1/2)`), or reuse `tol` for simplicity.

#### Constrained VV driver

```
# B
v  ← v + (dt/2) * a(r)
# A
r  ← r + dt * v
# SHAKE
apply_shake!(ps, cons, dt)
# recompute forces at corrected r
a  ← F(r) ./ m
# B
v  ← v + (dt/2) * a
# RATTLE (velocities)
apply_rattle!(ps, cons)
```

This driver preserves the geometric constraints to within tolerances and is compatible with the current force callback.

> **Thermostat compatibility:** If combining with Langevin (BAOAB), the O-step (stochastic OU) should be followed by a **RATTLE projection** of velocities to keep constraints satisfied. For this first iteration, we ship the **constrained VV** path; extending BAOAB with RATTLE is a follow-up.

### degrees\_of\_freedom with constraints

Let `N`, `D`, and `C = length(cons.r0)`. Then

```
dof = N*D - C
if remove_com
   dof -= D
end
dof = max(dof, 0)
```

This integrates with `instantaneous_temperature(ps; kB)`.

### remove\_com\_motion!

Mass-weighted COM velocity:

```
Vcom = (Σ m_i v_i) / (Σ m_i)
v_i ← v_i - Vcom
```

Optionally also recenter positions similarly (for non-PBC use).

---

## Numerical Pitfalls

* **Convergence:** Tight `tol` with complex constraint networks (e.g., rings) can require many iterations. Expose `maxiter` and return a clear error if not converged.
* **Ill-conditioning:** Very small `σ_l` (near-zero bond length or coincident atoms) leads to numerical blow-up; detect and abort early.
* **PBC displacement:** Minimum-image must be used **consistently** for `d_ij`. If constraints span a molecule that crosses a boundary, ensure unwrapped molecular coordinates are conceptually consistent. For now we assume constraints occur within a molecule and the minimum-image choice is valid.
* **Thermostats:** Any velocity randomization step must be followed by **RATTLE** to avoid drift off the constraint manifold.
* **DoF accounting:** When constraints are active (and optionally COM removal), temperature estimators and barostats must use the **reduced DoF** to avoid biased thermodynamics.
* **Large dt:** Constraints allow larger `dt` but do not make dynamics unconditionally stable; monitor energy and constraint residuals.

---

## Acceptance Tests

Add `test/test_constraints.jl` and include from `test/runtests.jl`.

```julia
@testset "DistanceConstraints basics" begin
    # Diatomic, 3D, one constraint
    N, D = 2, 3
    r0 = 1.25
    ps = ParticleSystem(zeros(N,D), zeros(N,D), ones(N))
    ps.positions[1,1] = 0.0
    ps.positions[2,1] = r0 + 0.1               # slight violation
    cons = DistanceConstraints([(1,2)], [r0]; tol=1e-10, maxiter=100)

    apply_shake!(ps, cons, 0.01)
    d = ps.positions[1,:] .- ps.positions[2,:]
    @test isapprox(norm(d), r0; rtol=0, atol=1e-8)
end)

@testset "RATTLE projects velocities" begin
    N, D = 2, 3
    r0 = 1.0
    ps = ParticleSystem(zeros(N,D), zeros(N,D), ones(N))
    ps.positions[1,1] = 0.0
    ps.positions[2,1] = r0
    cons = DistanceConstraints([(1,2)], [r0])

    # Give violating relative velocity along the bond
    ps.velocities[1,1] =  +0.3
    ps.velocities[2,1] =  -0.1
    apply_rattle!(ps, cons)
    d = ps.positions[1,:] .- ps.positions[2,:]
    vrel = ps.velocities[1,:] .- ps.velocities[2,:]
    @test isapprox(dot(d, vrel), 0.0; atol=1e-10)
end)

@testset "Constrained VV holds bond length" begin
    # Zero external forces: constraint should keep distance fixed over many steps
    N, D = 2, 3
    r0 = 0.75
    ps = ParticleSystem(zeros(N,D), zeros(N,D), ones(N))
    ps.positions[2,1] = r0
    ps.velocities .= 0.01                         # small motion
    cons = DistanceConstraints([(1,2)], [r0]; tol=1e-10, maxiter=100)
    dt = 0.02
    forces(R) = zero(R)                           # no forces

    for _ in 1:500
        velocity_verlet_shake_rattle!(ps, forces, dt, cons)
    end
    d = ps.positions[1,:] .- ps.positions[2,:]
    @test isapprox(norm(d), r0; atol=1e-7)
end)

@testset "degrees_of_freedom with constraints and COM" begin
    N, D = 5, 3
    pairs = [(1,2), (3,4)]
    r0s = [1.0, 1.5]
    cons = DistanceConstraints(pairs, r0s)
    ps = ParticleSystem(zeros(N,D), zeros(N,D), ones(N))
    @test degrees_of_freedom(ps; constraints=cons, remove_com=false) == N*D - length(pairs)
    @test degrees_of_freedom(ps; constraints=cons, remove_com=true)  == N*D - length(pairs) - D
end)

@testset "remove_com_motion! removes mass-weighted COM velocity" begin
    N, D = 3, 3
    ps = ParticleSystem(zeros(N,D), zeros(N,D), [1.0, 2.0, 3.0])
    ps.velocities .= 1.0
    remove_com_motion!(ps; which=:velocity)
    Vcom = (sum(ps.masses .* ps.velocities[:,1])) / sum(ps.masses)
    @test isapprox(Vcom, 0.0; atol=1e-14)
end)
```

---

## Implementation Notes

* Use `@views` and in-place updates; avoid temporary allocations in inner loops.
* Iterate constraints Gauss–Seidel style for better convergence; 2–5 sweeps usually suffice for simple molecules.
* Factor out a helper that computes `d_ij` (with minimum image if `use_minimum_image=true`) to keep code DRY.
* Keep `tol` on the squared constraint `|C_l|` (units length²). When checking `||r_i - r_j||`, convert appropriately.
* Extend `degrees_of_freedom` method to accept `constraints` and `remove_com` kwargs; keep old method for backward compatibility.

---

## Task for Implementer (small & focused)

1. Create `src/constraints.jl` with:

   * `struct DistanceConstraints` (fields as above) + constructor.
   * `apply_shake!`, `apply_rattle!` implementing the formulas and iterations described.
   * `velocity_verlet_shake_rattle!` that wraps existing VV steps and calls SHAKE/RATTLE appropriately.
   * `remove_com_motion!` (mass-weighted).
2. Extend `degrees_of_freedom(ps; constraints=nothing, remove_com=false)` in `src/thermostats.jl` (or a shared util) to account for constraints/COM.
3. `include("constraints.jl")` and export new symbols in `src/Verlet.jl`.
4. Add `test/test_constraints.jl` with the Acceptance Tests above and include from `test/runtests.jl`.
5. Ensure docstrings reference units and PBC options; add method signatures to the main README brief.

---

## To-Do’s for Documenter

* **Guide → Constraints:** Short page “Constrained Dynamics (SHAKE/RATTLE)” with:

  * When to use constraints; typical `dt` gains and caveats.
  * Minimal example constructing `DistanceConstraints` for a small molecule and running `velocity_verlet_shake_rattle!`.
  * Interplay with temperature: explain the **reduced DoF**.
* **API Reference:** `DistanceConstraints`, `apply_shake!`, `apply_rattle!`, `velocity_verlet_shake_rattle!`, `remove_com_motion!`, updated `degrees_of_freedom`.
* **Numerics Notes:** Tips on convergence (`tol`, `maxiter`), and PBC minimum image implications for intramolecular constraints.

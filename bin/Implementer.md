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
# DESIGN.md — Next Feature Plan

## Overview — **Constraint-Aware BAOAB (cBAOAB) + Residual Monitoring**

We now have working SHAKE/RATTLE and an unconstrained BAOAB thermostat. The natural next step is a **constrained Langevin integrator** that preserves holonomic distance constraints while sampling NVT:

* Implement **cBAOAB**: a BAOAB step that inserts **SHAKE** after position updates and **RATTLE** after velocity updates so that constraints remain satisfied despite stochastic kicks.
* Add a light utility for **constraint residual monitoring** (max norm and RMS) to aid debugging and regression tests.
* Keep scope tight: pairwise **DistanceConstraints** only (as already implemented), orthorhombic box support optional but consistent with minimum-image logic already present.

Benefits:

* Correct NVT sampling for systems with rigid bonds (e.g., water) while keeping the same user-facing force callback.
* Compatibility with existing `DistanceConstraints` and DoF accounting.

---

## Public API

```julia
"""
    langevin_baoab_constrained!(ps::ParticleSystem, forces, dt, cons::DistanceConstraints;
                                γ, T, kB=1.0, rng=Random.default_rng())

Advance one **constrained NVT** BAOAB step with SHAKE/RATTLE projections.
"""
langevin_baoab_constrained!(ps, forces, dt, cons; γ, T, kB=1.0, rng=Random.default_rng())

"""
    constraint_residuals(ps, cons::DistanceConstraints) -> (; maxC, rmsC, maxCd, rmsCd)

Return current position and velocity constraint residuals:
- C_l = ||r_i - r_j||^2 - r0_l^2
- Ċ_l = 2 (r_i - r_j) ⋅ (v_i - v_j)
"""
constraint_residuals(ps, cons)
```

**Exports (add):**

```julia
export langevin_baoab_constrained!, constraint_residuals
```

**Files (updated):**

* `src/thermostats.jl` — add `langevin_baoab_constrained!`.
* `src/constraints.jl` — add `constraint_residuals`.
* `src/Verlet.jl` — export new symbols.

---

## Data Structures

No new structs. Reuse:

```julia
DistanceConstraints  # immutable parameters; fields already implemented
ParticleSystem       # positions, velocities, masses
```

---

## Algorithms

### Constrained BAOAB (cBAOAB)

We follow the BAOAB splitting with **projections**:

Let `B` = half kick (deterministic), `A` = half drift (positions), `O` = OU stochastic velocity step.

We insert:

* **After each A:** `apply_shake!` (positions)
* **After each B or O:** `apply_rattle!` (velocities)

Step:

```
Inputs: ps = (R, V, m), forces, dt, γ, T, kB, rng, cons

1) B (half kick):           V ← V + (dt/2) * F(R) ./ m
   RATTLE (velocities):     apply_rattle!(ps, cons)          # keep Ċ = 0

2) A (half drift):          R ← R + (dt/2) * V
   SHAKE (positions):       apply_shake!(ps, cons, dt/2)

3) O (OU stochastic):
   c = exp(-γ*dt)
   for each particle i, component k:
       V[i,k] ← c * V[i,k] + sqrt((1-c^2) * kB*T / m_i) * ξ
   RATTLE (velocities):     apply_rattle!(ps, cons)

4) A (half drift):          R ← R + (dt/2) * V
   SHAKE (positions):       apply_shake!(ps, cons, dt/2)

5) Recompute forces:        F ← forces(R)

6) B (half kick):           V ← V + (dt/2) * F ./ m
   RATTLE (velocities):     apply_rattle!(ps, cons)

Return ps
```

Notes:

* We reuse the existing `_ou_coeff` small-x safeguard (`c ≈ 1 − γdt`) as in the unconstrained BAOAB.
* The `dt` argument passed to `apply_shake!` is not mathematically required but can be used for diagnostics; we keep the signature consistent with the existing function.
* Convergence tolerances and iteration caps come from `cons`.

### Residuals

For each constraint `l` with indices `(i,j)` and displacement `d_ij = r_i − r_j` (with minimum image if enabled):

```
C_l  = dot(d_ij, d_ij) − r0_l^2
Ċ_l = 2 * dot(d_ij, v_i − v_j)
```

Aggregate:

```
maxC  = maximum(abs, C_l)
rmsC  = sqrt(mean(C_l.^2))
maxCd = maximum(abs, Ċ_l)
rmsCd = sqrt(mean(Ċ_l.^2))
```

Return as a named tuple.

---

## Numerical Pitfalls

* **Projection frequency:** Omitting RATTLE after `O` will cause constraint drift in velocities. The proposed placements keep both `C` and `Ċ` small.
* **Tol/γ interaction:** Larger friction and temperature produce larger raw OU kicks; ensure `cons.tol` and `maxiter` are adequate to re-project. Users can relax `tol` slightly for high-T runs if needed.
* **Mass variance:** RATTLE uses `1/m_i` weights; extremely small masses (e.g., constrained hydrogens) can slow convergence if constraints are highly coupled. This is expected and acceptable for simple pair constraints.
* **DoF accounting:** `instantaneous_temperature` currently calls `degrees_of_freedom(ps)` (no kwargs). For constrained NVT runs, **document** that users should compute temperature using
  `degrees_of_freedom(ps; constraints=cons, remove_com=false)` if they want manual diagnostics. (We keep the API stable; no hidden global state.)
* **Random number reproducibility:** OU step uses `rng`. Tests should set a seed.

---

## Acceptance Tests

Add `test/test_cbaoab.jl` and include from `test/runtests.jl`.

```julia
@testset "cBAOAB preserves constraints (zero forces)" begin
    N, D = 3, 3
    r0 = 1.0
    ps = ParticleSystem(zeros(N,D), zeros(N,D), ones(N))
    ps.positions[2,1] = r0
    cons = DistanceConstraints([(1,2)], [r0]; tol=1e-10, maxiter=200)

    forces(R) = zero(R)
    dt, γ, T = 0.01, 1.0, 0.5
    rng = MersenneTwister(1234)

    # Warmup and integrate
    for _ in 1:1000
        langevin_baoab_constrained!(ps, forces, dt, cons; γ=γ, T=T, rng=rng)
    end

    d = ps.positions[1,:] .- ps.positions[2,:]
    @test isapprox(norm(d), r0; atol=1e-6)

    # Velocity residual should be tiny
    (; maxC, maxCd) = constraint_residuals(ps, cons)
    @test maxC ≤ 1e-8
    @test maxCd ≤ 1e-8
end)

@testset "cBAOAB temperature sane and finite" begin
    N, D = 8, 3
    ps = ParticleSystem(randn(N,D), zeros(N,D), ones(N))
    cons = DistanceConstraints([(1,2)], [1.0])
    forces(R) = zero(R)
    dt, γ, T = 0.005, 2.0, 1.2
    rng = MersenneTwister(7)

    # Run; measure running average temperature
    accT = 0.0
    steps = 5000
    for n in 1:steps
        langevin_baoab_constrained!(ps, forces, dt, cons; γ=γ, T=T, rng=rng)
        accT += instantaneous_temperature(ps; kB=1.0)
    end
    Tbar = accT / steps
    @test isfinite(Tbar)
    # Loose bound (small system; not asserting equality to T)
    @test 0.2 ≤ Tbar ≤ 3.0
end)

@testset "constraint_residuals reports zero at exact satisfaction" begin
    ps = ParticleSystem([0.0 0 0; 1.0 0 0], zeros(2,3), ones(2))
    cons = DistanceConstraints([(1,2)], [1.0])
    (; maxC, rmsC, maxCd, rmsCd) = constraint_residuals(ps, cons)
    @test maxC == 0.0
    @test rmsC == 0.0
    @test maxCd == 0.0
    @test rmsCd == 0.0
end)
```

---

## Implementation Notes

* Place the projections exactly as specified; avoid “batching” multiple OU kicks between projections.
* Reuse allocation-free patterns (`@views`, in-place ops). Use the same `invm = 1.0 ./ m` approach as in existing code.
* `constraint_residuals` must honor `use_minimum_image` and share the `_displacement!` helper from `constraints.jl` to avoid code duplication.

---

## Task for Implementer (small & focused)

1. **Add** `langevin_baoab_constrained!` to `src/thermostats.jl` following the algorithm above (B/half → RATTLE → A/half → SHAKE → O → RATTLE → A/half → SHAKE → recompute F → B/half → RATTLE).
2. **Add** `constraint_residuals(ps, cons)` to `src/constraints.jl` using the same minimum-image logic as `_displacement!`.
3. **Wire up** exports in `src/Verlet.jl`.
4. **Tests:** create `test/test_cbaoab.jl` with the Acceptance Tests and include it from `test/runtests.jl`.
5. **Docs:** briefly note in the BAOAB docstring that with constraints one should use `langevin_baoab_constrained!`.

---

## To-Do’s for Documenter

* **Guide → Constraints:** Add a short subsection “Constrained NVT (cBAOAB)” explaining where the projections are applied and why.
* **API Reference:** Document `langevin_baoab_constrained!` and `constraint_residuals`.
* **Numerics Notes:** Brief paragraph on the interplay of `γ`, `T`, and `cons.tol` and recommended tolerances.

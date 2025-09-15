# Spec: Module `Verlet.Constraints`

Purpose: Represent and enforce pair distance constraints via SHAKE/RATTLE, provide diagnostics and utilities.

## Types

- `struct DistanceConstraints{T_Int,T_Float}`
  - `i::Vector{T_Int}`, `j::Vector{T_Int}`: constrained pairs (1-based indices).
  - `r0::Vector{T_Float}`: target distances.
  - `tol::T_Float`, `maxiter::T_Int`: solver tolerance and max iterations.
  - `use_minimum_image::Bool`: whether to apply minimum-image using `sys.box`.
  - Convenience constructor: `DistanceConstraints(pairs::Vector{<:Tuple}, lengths::Vector{<:Real}; tol=1e-8, maxiter=50, use_minimum_image=true)`.

## Solvers

- `apply_shake!(sys::System, cons::DistanceConstraints, dt)`
  - Iteratively correct positions to satisfy `‖r_i − r_j‖^2 = r0^2` up to `tol`.
  - Uses Lagrange multiplier update `Δλ = −C/(2σ)` with `σ = (1/m_i + 1/m_j)‖Δ‖^2`.
  - Errors if not converged within `maxiter` or if `σ ≈ 0` (ill-conditioned).
- `apply_rattle!(sys::System, cons::DistanceConstraints)`
  - Correct velocities to satisfy `Ċ = 2Δ⋅(v_i − v_j) = 0` up to `tol`.
  - Similar iteration with `μ = −(Δ⋅v_rel) / τ`, `τ = (1/m_i + 1/m_j)‖Δ‖^2`.

## Time integration

- `velocity_verlet_shake_rattle!(sys, forces, dt, cons)`
  - Constrained velocity-Verlet: half kick → drift → SHAKE → new forces → half kick → RATTLE.

## Utilities

- `remove_com_motion!(sys; which=:velocity|:position|:both)`
  - Remove center-of-mass motion from velocities and/or positions.
- `constraint_residuals(sys::System, cons::DistanceConstraints)` → `(; maxC, rmsC, maxCd, rmsCd)`
  - Reports maximum and RMS residuals for positions and velocities.

## Performance

- Each SHAKE/RATTLE sweep is `O(#constraints)`; iteration count depends on tolerance and system state.
- `use_minimum_image` reduces discontinuities under PBCs for bonds crossing boundaries.

## Failure modes

- Ill-conditioned constraints (zero-length Δ or huge mass mismatch) may trigger errors.
- Non-convergence raises an error with the configured `maxiter`.

## Example

```julia
using Verlet, StaticArrays
R = [SVector(0.0,0,0), SVector(1.0,0,0)]
V = [SVector(0.0,0,0), SVector(0.0,0,0)]
F = [SVector(0.0,0,0), SVector(0.0,0,0)]
sys = System(R, V, F, [1.0,1.0], CubicBox(10.0), [1,1], Dict(1=>:A))
cons = DistanceConstraints([(1,2)], [1.0])
velocity_verlet_shake_rattle!(sys, R->F, 0.001, cons)
res = constraint_residuals(sys, cons)
```


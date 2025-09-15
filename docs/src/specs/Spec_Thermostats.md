# Spec: Module `Verlet.Thermostats`

Purpose: Thermostatting utilities and Langevin BAOAB integrators (un/constrained).

## Degrees of freedom and temperature

- `degrees_of_freedom(sys; constraints=nothing, remove_com=false) -> Int`
  - `N*D` reduced by number of constraints and optionally by `D` if COM removed.
- `instantaneous_temperature(sys; kB=1.0)`
  - `T = 2*KE / (kB * dof)` using `degrees_of_freedom(sys)`.
- `velocity_rescale!(sys, T; kB=1.0)`
  - Deterministically rescales velocities by `λ = sqrt(T/ max(Tinst, eps()))`.

## Langevin BAOAB (unconstrained)

- `langevin_baoab!(sys, forces, dt; γ, temp, kB=1.0, rng=Random.default_rng())`
  - Splitting: B(½) → A(½) → O → A(½) → B(½).
  - OU step uses `c = exp(-γ*dt)` (with a stable series fallback internally) and noise variance `(1-c^2) * kB*temp / m_i` per component.
  - With `γ=0` reduces to deterministic velocity-Verlet.

## Constrained Langevin BAOAB

- `langevin_baoab_constrained!(sys, forces, dt, cons::DistanceConstraints; γ, temp, kB=1.0, rng=Random.default_rng())`
  - Splitting with projections:
    1. B(½) → `apply_rattle!`
    2. A(½) → `apply_shake!`
    3. O → `apply_rattle!`
    4. A(½) → `apply_shake!`
    5. Recompute forces
    6. B(½) → `apply_rattle!`
  - Recommended path when rigid bonds/constraints are present.

## Invariants & Notes

- All updates are in-place on `sys.positions`/`sys.velocities`.
- `forces(R)` must return a vector of force vectors aligned with `R` and is called twice per BAOAB step (unconstrained) or once before/after projections (constrained).
- For reproducibility, pass an explicit RNG.

## Example

```julia
using Verlet, StaticArrays, Random
box = CubicBox(5.0)
R = [@SVector randn(3) for _ in 1:8]; wrap_positions!(R, box)
sys = System(R, fill(@SVector zeros(3),8), fill(@SVector zeros(3),8), ones(8), box, ones(Int,8), Dict(1=>:A))
forces = R -> [ -r for r in R ]
langevin_baoab!(sys, forces, 0.001; γ=1.0, temp=1.0, rng=MersenneTwister(1))
```


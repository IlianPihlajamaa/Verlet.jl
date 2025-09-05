# Constraints

Molecular systems often contain rigid bonds or fixed distances between atoms.
These are represented using the [`DistanceConstraints`](@ref) type and enforced
using SHAKE (positions) and RATTLE (velocities).

## Residual Monitoring

The helper function [`constraint_residuals`](@ref) reports how well constraints are
satisfied at a given state:

* `maxC`, `rmsC` — maximum and RMS positional residuals
* `maxCd`, `rmsCd` — maximum and RMS velocity residuals

```@example
ps = ParticleSystem([0.0 0 0; 1.0 0 0], zeros(2,3), ones(2))
cons = DistanceConstraints([(1,2)], [1.0])
constraint_residuals(ps, cons)
```

These values are useful for debugging and regression testing.

## Usage with cBAOAB

When running constrained Langevin dynamics with
[`langevin_baoab_constrained!`](@ref), residuals should remain close to machine
precision (typically `1e-8` or smaller) if the solver tolerance is sufficiently
strict.

---

See also:

* [`apply_shake!`](@ref)
* [`apply_rattle!`](@ref)
* [`langevin_baoab_constrained!`](@ref)

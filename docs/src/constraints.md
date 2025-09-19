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
using Verlet, StaticArrays
positions = [SVector{3}(0.0, 0.0, 0.0), SVector{3}(1.0, 0.0, 0.0)]
velocities = [SVector{3}(0.0, 0.0, 0.0) for _ in 1:2]
forces = [SVector{3}(0.0, 0.0, 0.0) for _ in 1:2]
masses = ones(2)
box = CubicBox(10.0)
types = [1, 1]
type_names = Dict(1 => :A)
sys = System(positions, velocities, forces, masses, box, types, type_names)
cons = DistanceConstraints([(1,2)], [1.0])
constraint_residuals(sys, cons)
```

These values are useful for debugging and regression testing.

## Usage with cBAOAB

When running constrained Langevin dynamics with
[`LangevinBAOABConstrained`](@ref Verlet.Integrators.LangevinBAOABConstrained),
residuals should remain close to machine precision (typically `1e-8` or
smaller) if the solver tolerance is sufficiently strict.

---

See also:

* [`apply_shake!`](@ref)
* [`apply_rattle!`](@ref)
* [`LangevinBAOABConstrained`](@ref Verlet.Integrators.LangevinBAOABConstrained)

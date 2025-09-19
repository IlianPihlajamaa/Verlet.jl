# Constrained Dynamics

Distance constraints lock specific inter-particle separations, enabling longer
timesteps and rigid fragments. Verlet.jl implements the classic SHAKE/RATTLE
family and integrates with both plain velocity Verlet and the Langevin
thermostat.

## Defining constraints

[`DistanceConstraints`](@ref Verlet.Constraints.DistanceConstraints) maps atom-index pairs to target lengths and solver tolerances.

Create a constraint set by listing atom index pairs and their target lengths:

```@example constraints_guide
using StaticArrays, LinearAlgebra, Verlet

pairs   = [(1, 2), (2, 3)]
lengths = [1.0, 1.0]
cons = DistanceConstraints(pairs, lengths; tol = 1e-10, maxiter = 100)
```

- `tol` controls the acceptable squared residual `|‖r_i - r_j‖² - r₀²|`.
- `maxiter` bounds the number of SHAKE/RATTLE correction sweeps per timestep.
- Set `use_minimum_image = true` when constraints may straddle periodic
  boundaries.

## Applying constraints during integration

Two options are available:

1. [`velocity_verlet_shake_rattle!`](@ref Verlet.Constraints.velocity_verlet_shake_rattle!)
   wraps the standard velocity-Verlet update with SHAKE (position projection) and
   RATTLE (velocity projection).
2. [`Verlet.Integrators.LangevinBAOABConstrained`](@ref) embeds the same
   projections in the stochastic BAOAB scheme.

Both require a callable that returns forces for the updated positions. The force
function should honour the `System` stored in the closure.

```@example constraints_guide
using StaticArrays, LinearAlgebra, Verlet

N = 2
zero3 = SVector{3}(0.0, 0.0, 0.0)
positions = [zero3, SVector{3}(1.05, 0.0, 0.0)]
velocities = fill(zero3, N)
forces = fill(zero3, N)
masses = ones(N)
box = CubicBox(10.0)
types = ones(Int, N)
type_names = Dict(1 => :A)
ff = ForceField(())  # empty forcefield placeholder
sys = System(positions, velocities, forces, masses, box, types, type_names;
             forcefield = ff)

forces_fn(R) = begin
    _ = R
    Verlet.Core.compute_all_forces!(sys)
    sys.forces
end

dt = 0.01
velocity_verlet_shake_rattle!(sys, forces_fn, dt, DistanceConstraints([(1, 2)], [1.0]))
round(norm(sys.positions[1] - sys.positions[2]), digits = 6)
```

## Diagnostics and utilities

- [`constraint_residuals`](@ref Verlet.Constraints.constraint_residuals)
  returns maximum and RMS residuals for both positions and velocities.
- [`remove_com_motion!`](@ref Verlet.Constraints.remove_com_motion!) eliminates
  centre-of-mass drift in velocity or position space, useful prior to applying
  thermostats or barostats.
- [`Verlet.Thermostats.degrees_of_freedom`](@ref) accepts the `constraints`
  keyword so temperature estimators account for removed degrees of freedom.

## Choosing tolerances

- `tol = 1e-8` suffices for most biomolecular applications; tighten to `1e-10`
  for high-precision energy conservation.
- Increase `maxiter` when solving strongly coupled constraint networks (e.g.
  rings). Failure to converge raises an error so you can respond appropriately.

## Interplay with thermostats

Thermostat updates that randomise velocities (e.g. BAOAB) must immediately be
followed by RATTLE projection. The constrained Langevin integrator handles this
internally; for custom schemes, call `apply_rattle!(sys, cons)` after modifying
velocities.

## Further reading

- Ryckaert, Ciccotti & Berendsen (1977). *Numerical integration of the Cartesian
  equations of motion of a system with constraints.*
- Andersen (1983). *RATTLE: A "velocity" version of the SHAKE algorithm.*

For a step-by-step build, revisit
[Tutorial 3 · Constraints in Practice](../tutorials/constraints.md).

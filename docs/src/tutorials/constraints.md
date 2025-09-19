# Tutorial 3 · Constraints in Practice

Rigid constraints let you hold specific distances fixed, enabling larger
integration timesteps and faithful rigid fragments (e.g. water molecules). This
tutorial uses Verlet.jl's SHAKE/RATTLE implementation to maintain a fixed bond.

## Scenario

Simulate a dumbbell molecule with a 1.0 distance between two atoms, using a
custom force and the `velocity_verlet_shake_rattle!` driver.

## 1. Define the model

```@example constrained
using LinearAlgebra, StaticArrays, Verlet

struct Tether
    k::Float64
end

function Verlet.Core.compute_forces!(pot::Tether, sys::System)
    # Simple harmonic tether to the origin for each atom
    @inbounds for i in 1:natoms(sys)
        sys.forces[i] += -pot.k * sys.positions[i]
    end
    return sys
end
```

## 2. Create the constrained system

```@example constrained
zero3 = SVector{3}(0.0, 0.0, 0.0)
positions = [zero3, SVector{3}(1.1, 0.0, 0.0)]
velocities = [zero3, zero3]
forces = [zero3, zero3]
masses = [1.0, 1.0]
box = CubicBox(10.0)
types = [1, 1]
type_names = Dict(1 => :A)
cons = DistanceConstraints([(1, 2)], [1.0]; tol = 1e-10, maxiter = 100)
ff = ForceField((Tether(0.2),))
sys = System(positions, velocities, forces, masses, box, types, type_names;
             forcefield = ff, specific_potentials = ())
```

## 3. Advance with SHAKE/RATTLE

The helper `velocity_verlet_shake_rattle!` performs a velocity-Verlet step and
projects onto the constraint manifold each sub-step.

```@example constrained
using LinearAlgebra

forces_fn(R) = begin
    _ = R  # system.positions already updated in-place
    Verlet.Core.compute_all_forces!(sys)
    sys.forces
end

dt = 0.01
for step in 1:200
    velocity_verlet_shake_rattle!(sys, forces_fn, dt, cons)
end

round(norm(sys.positions[1] - sys.positions[2]), digits = 6)
```

The distance is now maintained at 1.0 (within tolerance) throughout the
trajectory.

## 4. Monitoring residuals

`constraint_residuals` provides diagnostics for both positional and velocity
violations—essential for validating convergence when you relax tolerances or
increase the timestep.

```@example constrained
Verlet.Constraints.constraint_residuals(sys, cons)
```

## Checklist for stable constrained runs

- Select a tight enough tolerance for your system; `1e-8–1e-10` is typical.
- Increase `maxiter` if the residuals fail to drop below `tol`.
- Pass the same `cons` object to any thermostats so velocities remain
  constraint-consistent after stochastic steps.
- Use `degrees_of_freedom(sys; constraints=cons, remove_com=true)` when
  computing temperatures.

Next, continue with [Tutorial 4 · Thermostatted Dynamics](thermostat.md) to add a
Langevin thermostat on top of constrained motion.

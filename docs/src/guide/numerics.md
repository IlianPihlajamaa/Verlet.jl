# Numerical Notes & Best Practices

This page collects practical advice for choosing timesteps, monitoring energy,
and combining integrators with thermostats in Verlet.jl.

## Timestep heuristics

- Start by scanning the largest timestep that conserves total energy for an
  unconstrained, microcanonical (NVE) run. Monitor `kinetic_energy + potential`
  over several hundred steps.
- Holonomic constraints (SHAKE/RATTLE) allow roughly 2× larger timesteps compared
  with fully flexible bonds.
- Stiff bonded terms or strong Lennard-Jones overlaps may necessitate smaller
  timesteps; consider energy minimisation (e.g. conjugate gradient) before
  production dynamics.

## Energy monitoring pattern

```@example numerics
using StaticArrays, Verlet

struct Hooke
    k::Float64
end

function Verlet.Core.compute_forces!(pot::Hooke, sys::System)
    @inbounds for i in 1:natoms(sys)
        sys.forces[i] += -pot.k * sys.positions[i]
    end
    return sys
end

box = CubicBox(10.0)
positions = [SVector{3}(1.0, 0.0, 0.0)]
velocities = [SVector{3}(0.0, 1.0, 0.0)]
forces = [SVector{3}(0.0, 0.0, 0.0)]
masses = [1.0]
types = [1]
type_names = Dict(1 => :A)
hooke = Hooke(1.0)
ff = ForceField((hooke,))
sys = System(positions, velocities, forces, masses, box, types, type_names;
             forcefield = ff)

vv = VelocityVerlet(0.05)
energies = Float64[]
integrate!(vv, sys, 400; callback = (system, step, _) -> begin
    U = 0.5 * hooke.k * sum(abs2, system.positions[1])
    push!(energies, kinetic_energy(system) + U)
    return nothing
end)

(round(minimum(energies), digits = 6), round(maximum(energies), digits = 6))
```

Use a similar callback to validate new potentials or parameter sets.

## Thermostats

- `Verlet.Integrators.LangevinBAOAB` is a good default for canonical (NVT)
  sampling; it preserves configurational accuracy and requires only one force
  evaluation per step.
- Set `γ*dt` between `0.1` and `1.0`. Smaller values retain more ballistic
  character; larger values overdamp dynamics but still sample the Boltzmann
  distribution.
- Always supply an explicit RNG (`rng = MersenneTwister(seed)`) for reproducible
dynamics.
- For constrained systems, switch to `LangevinBAOABConstrained` with the same
  `DistanceConstraints` object used by SHAKE/RATTLE.

## Temperature estimation

```
T_inst = 2 * kinetic_energy(sys) / (kB * dof)
```

Use `Verlet.Thermostats.instantaneous_temperature(sys; kB)` for convenience.
Remember to adjust the degrees of freedom:

```@example numerics
cons = DistanceConstraints([(1,2)], [1.0])
Verlet.Thermostats.degrees_of_freedom(sys; constraints = cons, remove_com = true)
```

## Random initial velocities

Draw initial velocities compatible with your target temperature:

```@example numerics
using Random
Random.seed!(123)
Ttarget = 1.5
for i in eachindex(sys.velocities)
    sys.velocities[i] = sqrt(Ttarget) * SVector{3}(randn(), randn(), randn())
end
Verlet.Thermostats.remove_com_motion!(sys; which = :velocity)
```

## Troubleshooting

- **Exploding energies:** reduce the timestep, minimise initial structures, or
  double-check force signs. Large overlaps can produce huge forces.
- **Constraint failures:** diagnose with `constraint_residuals`. Raise
  `maxiter`, loosen `tol`, or reduce the timestep.
- **Neighbour rebuild storms:** increase the pair potential `skin` (e.g. from
  `0.3` to `0.6`) so particles travel farther before triggering a rebuild.

For a concrete thermostat example, see
[Tutorial 4 · Thermostatted Dynamics](../tutorials/thermostat.md).

## Integrator API

```@docs
Verlet.Integrators.VelocityVerlet
Verlet.Integrators.ConjugateGradient
Verlet.Integrators.LangevinBAOAB
Verlet.Integrators.LangevinBAOABConstrained
```

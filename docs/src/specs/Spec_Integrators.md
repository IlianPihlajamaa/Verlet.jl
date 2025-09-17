# Spec: Module `Verlet.Integrators`

Purpose: Provide concrete `AbstractIntegrator` implementations for common
molecular dynamics workflows.

## Types

```@docs
Verlet.Integrators.VelocityVerlet
Verlet.Integrators.ConjugateGradient
```

- `VelocityVerlet(dt; wrap=false)`
  - Fields: timestep `dt`, periodic wrapping flag.
  - `integrate!(VelocityVerlet, system, nsteps; callback=nothing)` performs
    velocity-Verlet integration using the `ForceField` attached to `system`;
    the callback receives `(system, step)` and may return `false` to halt early.
- `ConjugateGradient(energy; tol=1e-8, alpha0=1.0, min_alpha=1e-8, c1=1e-4, wrap=true)`
  - Implements non-linear Polak–Ribière conjugate-gradient minimisation with
    Armijo backtracking. The `energy` callable must accept the `system` and
    return the scalar potential energy. `integrate!` treats `nsteps` as the
    maximum number of line-search iterations.

## Behaviour & Invariants

- `integrate!` on each integrator mutates the supplied `System` in place and
  always updates `system.forces` with the most recent evaluation.
- Callbacks returning `false` stop integration early without error.
- Negative `nsteps` raise `ArgumentError`.
- Periodic wrapping is optional (`wrap=true` wraps positions after each update).

## Example

```julia
using Verlet, StaticArrays

struct Springs
    k::Float64
end

function Verlet.Core.compute_forces!(pot::Springs, sys::System)
    @inbounds for i in 1:natoms(sys)
        sys.forces[i] += -pot.k * sys.positions[i]
    end
    return sys
end

box = CubicBox(10.0)
R = [@SVector randn(3) for _ in 1:4]; wrap_positions!(R, box)
sys = System(R, fill(@SVector zeros(3), 4), fill(@SVector zeros(3), 4), ones(4), box, ones(Int,4), Dict(1=>:A);
           forcefield=ForceField((Springs(1.0),)))
vv = VelocityVerlet(0.001)
integrate!(vv, sys, 100)
```

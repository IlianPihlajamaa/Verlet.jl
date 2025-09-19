# Spec: Module `Verlet.Thermostats`

Purpose: Thermostatting utilities (degrees of freedom, temperature estimators, deterministic rescaling).

## Degrees of freedom and temperature

- `degrees_of_freedom(sys; constraints=nothing, remove_com=false) -> Int`
  - `N*D` reduced by number of constraints and optionally by `D` if COM removed.
- `instantaneous_temperature(sys; kB=1.0)`
  - `T = 2*KE / (kB * dof)` using `degrees_of_freedom(sys)`.
- `velocity_rescale!(sys, T; kB=1.0)`
  - Deterministically rescales velocities by `λ = sqrt(T / max(Tinst, eps()))`.

## Langevin Integrators

Langevin BAOAB schemes now live in [`Verlet.Integrators`](Spec_Integrators.md):

- `LangevinBAOAB(dt; γ, temp, kB=1.0, wrap=false, rng=Random.default_rng())`
- `LangevinBAOABConstrained(dt, constraints; γ, temp, kB=1.0, wrap=false, rng=Random.default_rng())`

Each implements `step!(integrator, system)` and works with `integrate!`, rebuilding
neighbor lists and updating forces via the `System`'s `ForceField`.

## Example

```julia
using Verlet, StaticArrays, Random

struct Springs
    k::Float64
end

function Verlet.Core.compute_forces!(pot::Springs, sys::System)
    @inbounds for i in 1:natoms(sys)
        sys.forces[i] += -pot.k * sys.positions[i]
    end
    return sys
end

box = CubicBox(5.0)
R = [@SVector randn(3) for _ in 1:8]; wrap_positions!(R, box)
sys = System(R, fill(@SVector zeros(3),8), fill(@SVector zeros(3),8), ones(8), box, ones(Int,8), Dict(1=>:A);
           forcefield=ForceField((Springs(1.0),)))

integrator = Verlet.Integrators.LangevinBAOAB(0.001; γ=1.0, temp=1.0, rng=MersenneTwister(1))
integrate!(integrator, sys, 100)
```

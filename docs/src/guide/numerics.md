# NVT (Langevin BAOAB) Thermostat

This guide shows how to run **constant-temperature (NVT)** dynamics using the
`Verlet.Integrators.LangevinBAOAB` integrator and how to monitor/adjust the
system temperature with a few convenience tools.

## Why BAOAB?

- Excellent configurational sampling compared to naive Euler–Maruyama.
- Robust for larger `dt` at a given target `T`.
- Exact Ornstein–Uhlenbeck (OU) substep and only **one force evaluation per step**.
- With `γ = 0`, BAOAB reduces to velocity–Verlet.

## Quick Start

```julia
using Random, StaticArrays, Verlet

struct Springs
    k::Float64
end

function Verlet.Core.compute_forces!(pot::Springs, sys::System)
    @inbounds for i in 1:natoms(sys)
        sys.forces[i] += -pot.k * sys.positions[i]
    end
    return sys
end

N, D = 64, 3
positions = [@SVector zeros(D) for _ in 1:N]
velocities = [@SVector randn(D) for _ in 1:N]
forces = [@SVector zeros(D) for _ in 1:N]
masses = ones(N)
box = CubicBox(10.0)
types = ones(Int, N)
type_names = Dict(1 => :A)
ff = ForceField((Springs(0.2),))
ps = System(positions, velocities, forces, masses, box, types, type_names; forcefield=ff)

dt = 0.005      # time step
γ  = 1.0        # friction [1/time]
Tt = 1.5        # target temperature (kB = 1.0 here)
rng = MersenneTwister(2025)

# Optional: bring KE near the target quickly
Verlet.Thermostats.velocity_rescale!(ps, Tt; kB=1.0)

integrator = Verlet.Integrators.LangevinBAOAB(dt; γ=γ, temp=Tt, rng=rng)
integrate!(integrator, ps, 2_000)

@info "Instantaneous T" Verlet.Thermostats.instantaneous_temperature(ps; kB=1.0)
```

## Choosing `γ` and `dt`

- **Units:** ensure `γ` is in `1/time` and `kB*T` in energy, consistent with `r`, `v`, and `m`.
- Start with `γ*dt ≈ 0.1–1`. Larger values overdamp dynamics (OK for sampling; poor for kinetics).
- Keep `dt` similar to your velocity-Verlet choice at the same force field. BAOAB is often **at least as stable**.

## Algorithm Sketch

BAOAB splits each step into **B (half kick) → A (half drift) → O (OU) → A → B**. The OU substep uses the exact update
`v ← c*v + s*ξ` with `c = exp(-γ*dt)` and `s = √((1 - c^2) * kB*T / m)`. See the integrator docstring for further details,
including numerics for small `γ*dt`.

## Diagnostics

- `Verlet.Thermostats.instantaneous_temperature(ps; kB)` returns the current kinetic temperature using equipartition (`dof = N*D` for now).
- Time-average `Verlet.Thermostats.instantaneous_temperature` over many steps (or save to disk) to validate the thermostat.

## Pitfalls & Tips

- **Mass heterogeneity:** compute OU noise per particle using its mass; do **not** reuse a single scale across species.
- **Small systems:** temperature fluctuates strongly; judge convergence statistically.
- **Reproducibility:** always pass an explicit RNG when constructing `LangevinBAOAB`.
- **Not a thermostat:** `Verlet.Thermostats.velocity_rescale!` is useful for a single "KE nudge" but does not produce the correct canonical distribution.

## API Reference

```@docs
Verlet.Integrators.LangevinBAOAB
Verlet.Integrators.LangevinBAOABConstrained
Verlet.Thermostats.instantaneous_temperature
Verlet.Thermostats.degrees_of_freedom
Verlet.Thermostats.velocity_rescale!
```

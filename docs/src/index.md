# Verlet.jl

A lightweight toolkit for molecular dynamics experiments in Julia.
Verlet.jl offers composable force fields, neighbour-list acceleration, and
convenient drivers for both deterministic and thermostatted dynamics.

## Installation

From the Julia package REPL:

```julia
pkg> add Verlet
```

Verify the install and load the main module:

```julia
julia> using Verlet
```

## A ten-line simulation

```@example intro
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
positions  = [SVector{3}(1.0, 0.0, 0.0)]
velocities = [SVector{3}(0.0, 0.5, 0.0)]
forces     = [SVector{3}(0.0, 0.0, 0.0)]
masses     = [1.0]
types      = [1]
type_names = Dict(1 => :A)
ff = ForceField((Hooke(1.0),))
sys = System(positions, velocities, forces, masses, box, types, type_names; forcefield = ff)

vv = VelocityVerlet(0.05)
integrate!(vv, sys, 200)
sys.positions[1]
```

## Documentation map

- **Tutorials** – follow the numbered lessons: start with
  [Tutorial 1 · Hello, Velocity Verlet](tutorials/quickstart.md) and work your
  way towards pair potentials, constraints, and thermostats.
- **How-to Guides** – focused references for day-to-day tasks such as
  [building systems](guide/system.md),
  [defining forces](guide/forces.md), or
  [managing neighbour lists](guide/neighbors.md).
- **Reference** – specification sheets for each submodule and the generated
  [API reference](api.md) when you need precise type signatures.

## Highlights

- 🔁 Velocity Verlet, conjugate gradient, and Langevin BAOAB integrators.
- ⚛️ Built-in Lennard-Jones, Coulomb, and bonded potentials with exclusion
  support.
- 📦 Holonomic constraints (SHAKE/RATTLE) and thermostat-aware variants.
- 📚 First-class `StaticArrays` support for allocation-free inner loops.
- 🧮 Observable/logging framework with built-in thermodynamic probes.

Ready to dive in? Start with
[Building Systems](guide/system.md) or jump straight to the
[tutorial series](tutorials/quickstart.md).

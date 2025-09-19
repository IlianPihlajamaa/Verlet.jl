# Tutorial 1 · Hello, Velocity Verlet

This tutorial walks through the minimum ingredients required to step a tiny
molecular system with Verlet.jl. We will integrate a single particle tethered
by a harmonic spring and peek at the resulting trajectory and energies.

## Prerequisites

```julia
using Verlet, StaticArrays
```

The `StaticArrays` package is convenient for creating the `SVector`
coordinates used throughout Verlet.jl. All other functionality lives in the
`Verlet` module itself.

## 1. Define a force law

We implement a tiny Hookean force as a Julia struct plus a
`compute_forces!` overload. This mirrors how you would bring your own
potential to the framework.

```@example hello
using Verlet, StaticArrays

struct Hooke
    k::Float64
end

function Verlet.Core.compute_forces!(pot::Hooke, sys::System)
    @inbounds for i in 1:natoms(sys)
        sys.forces[i] += -pot.k * sys.positions[i]
    end
    return sys
end
```

## 2. Build a `System`

We supply positions, velocities, forces, masses, a periodic box, and particle
types. Everything is held in small static vectors for performance.

```@example hello
box = CubicBox(10.0)
zero3 = SVector{3}(0.0, 0.0, 0.0)
positions  = [SVector{3}(1.0, 0.0, 0.0)]
velocities = [SVector{3}(0.0, 0.5, 0.0)]
forces     = [zero3]
masses     = [1.0]
types      = [1]
type_names = Dict(1 => :A)

hooke = Hooke(1.0)
ff  = ForceField((hooke,))
sys = System(positions, velocities, forces, masses, box, types, type_names;
             forcefield = ff)
```

At this point we have a fully-specified mechanical system: positions and
velocities, an empty force accumulator, masses, and a force field containing
our Hookean potential.

## 3. Advance with velocity Verlet

Create an integrator, run for a few hundred steps, and capture the trajectory.

```@example hello
vv = VelocityVerlet(0.05)
trajectory = Vector{SVector{3,Float64}}()
integrate!(vv, sys, 400; callback = (system, step, _) -> begin
    push!(trajectory, system.positions[1])
    return nothing
end)

length(trajectory)
```

We can look at the final position and total energy:

```@example hello
final_r = trajectory[end]
final_E = kinetic_energy(sys) + 0.5 * hooke.k * sum(abs2, sys.positions[1])
(final_r, round(final_E, digits=6))
```

Try experimenting with the timestep, mass, or spring constant to see how the
motion changes.

## Where to go next

- [Tutorial 2 · Pair Potentials & Neighbour Lists](pair_potentials.md) introduces
  Lennard-Jones forces and automatic neighbour upkeep.
- Browse the [Building Systems](../guide/system.md) guide for a deeper tour of
the `System` container and associated helpers.

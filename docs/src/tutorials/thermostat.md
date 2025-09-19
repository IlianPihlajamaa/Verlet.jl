# Tutorial 4 · Thermostatted Dynamics

Langevin BAOAB integrates Newton's equations while gently steering the kinetic
energy toward a chosen temperature. In this tutorial we add a thermostat to the
Lennard-Jones system from Tutorial 2 and track the instantaneous temperature.

## 1. Reuse the Lennard-Jones setup

```@example thermo
using Random, StaticArrays, Statistics, Verlet

Random.seed!(2025)
N = 64
box = CubicBox(8.0)
randv() = SVector{3}(randn(), randn(), randn())
positions = [randv() for _ in 1:N]
wrap_positions!(positions, box)
velocities = [randv() for _ in 1:N]
zero3 = SVector{3}(0.0, 0.0, 0.0)
forces = fill(zero3, N)
masses = fill(1.0, N)
types = ones(Int, N)
type_names = Dict(1 => :A)

params = PairTable(fill(LJPair(1.0, 1.0, 2.5), (1, 1)))
lj = LennardJones(params, Tuple{T_Int,T_Int}[], 0.4)
ff = ForceField((lj,))
sys = System(positions, velocities, forces, masses, box, types, type_names;
             forcefield = ff)
```

## 2. Prepare the thermostat

```@example thermo
dt = 0.005
γ = 1.0
Ttarget = 1.5
rng = MersenneTwister(42)

integrator = Verlet.Integrators.LangevinBAOAB(dt; γ = γ, temp = Ttarget, rng = rng)
```

Optional: rescale velocities once to start near the target temperature.

```@example thermo
Verlet.Thermostats.velocity_rescale!(sys, Ttarget)
```

## 3. Integrate and record the temperature trace

```@example thermo
samples = Float64[]
integrate!(integrator, sys, 2_000; callback = (system, step, _) -> begin
    push!(samples, Verlet.Thermostats.instantaneous_temperature(system))
    return nothing
end)

(mean(samples), std(samples))
```

The temperature fluctuates around the target, as expected for Langevin
thermostats.

## 4. Constrained variant

If you are also running with holonomic constraints, use
`LangevinBAOABConstrained` and provide the same `DistanceConstraints` instance
used by your constraint solver.

```julia
cons = DistanceConstraints([(1, 2)], [1.0])
thermo_cons = Verlet.Integrators.LangevinBAOABConstrained(dt, cons;
    γ = γ, temp = Ttarget, rng = rng)
```

## Tips

- `γ*dt` controls the relaxation strength. Values in `[0.1, 1.0]` are typical.
- Pass an explicit `rng` to make runs reproducible.
- Always compute instantaneous temperatures with the appropriate number of
  degrees of freedom: `degrees_of_freedom(sys; constraints = cons, remove_com = true)`.

Congratulations—you now have a complete NVT workflow. Dive into the
[How-to Guides](../guide/system.md) for deeper control over system setup and
the available force-field components.

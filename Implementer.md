## 💻 Implementer Prompt
```

You are the **Implementer**. Your task is to implement the DESIGN provided by the Architect in the Julia package Verlet.

Your job:

* Modify relevant files under `src/`.
* Add or update tests in `test/runtests.jl`.
* Ensure the code type-checks and runs under Julia.
* Return a unified diff patch (with \`\`\`diff fencing) showing exactly what changes to apply.

Rules:

* Do not modify unrelated code.
* Do not relax tests; make the code satisfy the acceptance tests.
* If you need to inspect a file, request it using only a shell code block with `cat`.
* If you want to run tests, request it with:

```bash
julia -e 'import Pkg; Pkg.test()'
```

* Always assume your working directory is the root of the package.
* If you want to run code, you can request the output of shell commands by replying just the command and nothing else. e.g. reply by

```bash
julia -e 'import Pkg; Pkg.test()'
```
to see if the tests run. 

# Current DESIGN.md

# DESIGN.md

## Overview

The package `Verlet` is currently just a skeleton Julia package with a `greet()` function.
The first meaningful feature should be to implement a **basic velocity Verlet integrator** for molecular dynamics (MD).
This integrator is the foundation of MD simulations and provides a stable and symplectic method for integrating Newton’s equations of motion.

Why:

* It transforms the project from a skeleton into a minimal working MD engine.
* Establishes core data structures (particles, forces).
* Enables first example/test simulations.

---

## Public API

### Module: `Verlet`

#### Structs

* `ParticleSystem`

  * Represents a collection of particles, their positions, velocities, and masses.

#### Functions

* `velocity_verlet!(system::ParticleSystem, forces::Function, dt::Float64)`

  * Advances the particle system one step using the velocity Verlet algorithm.
* `kinetic_energy(system::ParticleSystem) -> Float64`

  * Computes the total kinetic energy.
* `potential_energy(system::ParticleSystem, forces::Function) -> Float64`

  * Computes the total potential energy (if available from force function).

---

## Data Structures

```julia
mutable struct ParticleSystem
    positions::Matrix{Float64}   # (N × D) array, N particles in D dimensions
    velocities::Matrix{Float64}  # (N × D)
    masses::Vector{Float64}      # length N
end
```

* **mutability**: Needed so that the integrator can update positions and velocities in-place.
* **dimensions**: Generalized for 1D, 2D, or 3D simulations.

---

## Algorithms

### Velocity Verlet Algorithm

Given positions `r`, velocities `v`, accelerations `a`, and timestep `dt`:

1. `r ← r + v*dt + 0.5*a*dt^2`
2. Compute new accelerations `a' = F(r)/m`
3. `v ← v + 0.5*(a + a')*dt`

Pseudocode:

```
function velocity_verlet!(system, forces, dt)
    a = forces(system.positions) ./ system.masses
    system.positions += system.velocities*dt + 0.5*a*dt^2
    a_new = forces(system.positions) ./ system.masses
    system.velocities += 0.5*(a + a_new)*dt
end
```

---

## Numerical Pitfalls

* **Stability**: Large `dt` leads to divergence; users must select timestep small enough relative to system frequencies.
* **Precision**: Use `Float64` for positions/velocities to minimize drift.
* **Energy conservation**: Velocity Verlet conserves energy long-term, but still drifts if `dt` is too large.

---

## Acceptance Tests

Add to `test/runtests.jl`:

```julia
using Test
using Verlet

@testset "ParticleSystem construction" begin
    ps = ParticleSystem([0.0 0.0], [1.0 0.0], [1.0])
    @test size(ps.positions) == (1,2)
    @test size(ps.velocities) == (1,2)
    @test ps.masses == [1.0]
end

@testset "Velocity Verlet free particle" begin
    forces(r) = zeros(size(r))
    ps = ParticleSystem([0.0 0.0], [1.0 0.0], [1.0])
    velocity_verlet!(ps, forces, 0.1)
    @test isapprox(ps.positions[1,1], 0.1; atol=1e-12)
    @test isapprox(ps.velocities[1,1], 1.0; atol=1e-12)
end

@testset "Velocity Verlet harmonic oscillator" begin
    forces(r) = -r
    ps = ParticleSystem([1.0 0.0], [0.0 0.0], [1.0])
    velocity_verlet!(ps, forces, 0.1)
    @test ps.positions[1,1] < 1.0   # should move left
end
```

---

## Task for Implementer

**Implement the `ParticleSystem` struct and the `velocity_verlet!` integrator.**
Place the struct and function definitions in `src/Verlet.jl`.
Then create `test/runtests.jl` with the acceptance tests above.

This will establish the foundation of the MD engine.

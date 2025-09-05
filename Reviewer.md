## 🔍 Reviewer Prompt
```

You are the **Reviewer**. Your job is to review patches proposed by the Implementer who implemented
features given by Design.md.

Your job:

* Check correctness of Julia code (types, multiple dispatch, array broadcasting).
* Verify algorithm correctness against the DESIGN.
* Identify numerical pitfalls (round-off errors, energy drift, unstable integration).
* Check that all acceptance tests in the DESIGN are covered.
* Respond with ACCEPT (if okay) or REJECT (with specific reasons + fixes).

Rules:

* If you want to see the state of a file, request it with a shell command using `cat`.
* If you want to run code, you can request the output of shell commands by replying just the command and nothing else. e.g. reply by

```bash
julia -e 'import Pkg; Pkg.test()'
```
to see if the tests run. 

* Always assume your working directory is the root of the package.

# Design 

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


# Implementation



*** Begin Patch
*** Update File: src/Verlet.jl
@@
 module Verlet
 
-greet() = print("Hello World!")
+export ParticleSystem, velocity_verlet!, kinetic_energy, potential_energy
+
+"""
+    ParticleSystem
+
+Represents a collection of particles, their positions, velocities, and masses.
+
+Fields:
+- `positions::Matrix{Float64}`: (N × D) array of positions
+- `velocities::Matrix{Float64}`: (N × D) array of velocities
+- `masses::Vector{Float64}`: length-N vector of masses
+"""
+mutable struct ParticleSystem
+    positions::Matrix{Float64}   # (N × D)
+    velocities::Matrix{Float64}  # (N × D)
+    masses::Vector{Float64}      # (N)
+end
+
+"""
+    velocity_verlet!(system::ParticleSystem, forces::Function, dt::Float64)
+
+Advance the system one step using the velocity Verlet integrator.
+
+The `forces` function must accept a positions matrix `(N × D)` and return
+an equally-sized matrix of forces.
+"""
+function velocity_verlet!(system::ParticleSystem, forces::Function, dt::Float64)
+    # Current accelerations a = F/m, broadcasting over rows (per particle)
+    a = forces(system.positions) ./ system.masses
+    # Update positions
+    system.positions .+= system.velocities .* dt .+ 0.5 .* a .* (dt^2)
+    # New accelerations from updated positions
+    a_new = forces(system.positions) ./ system.masses
+    # Update velocities
+    system.velocities .+= 0.5 .* (a .+ a_new) .* dt
+    return system
+end
+
+"""
+    kinetic_energy(system::ParticleSystem) -> Float64
+
+Compute the total kinetic energy `sum(½ m v²)` of the system.
+"""
+function kinetic_energy(system::ParticleSystem)::Float64
+    v2_per_particle = sum(abs2, eachrow(system.velocities)) # length N
+    return 0.5 * sum(system.masses .* v2_per_particle)
+end
+
+"""
+    potential_energy(system::ParticleSystem, forces::Function) -> Float64
+
+Attempt to compute total potential energy. This is only possible if the
+`forces` callable supports one of the following conventions:
+
+1) Returns a tuple `(F, U)` when called with positions, where `F` is the
+   force matrix and `U` is the total potential energy (Float64).
+2) Supports a keyword `return_potential=true` and returns `(F, U)`.
+
+If neither convention is supported, this function throws an error.
+"""
+function potential_energy(system::ParticleSystem, forces::Function)::Float64
+    # Try convention 1: plain call returns (F, U)
+    F_U = try
+        forces(system.positions)
+    catch
+        nothing
+    end
+    if F_U !== nothing && F_U isa Tuple && length(F_U) == 2
+        _, U = F_U
+        return float(U)
+    end
+    # Try convention 2: keyword to request potential
+    F_U_kw = try
+        forces(system.positions; return_potential=true)
+    catch
+        nothing
+    end
+    if F_U_kw !== nothing && F_U_kw isa Tuple && length(F_U_kw) == 2
+        _, U = F_U_kw
+        return float(U)
+    end
+    error("Force function does not provide potential energy; cannot compute potential_energy.")
+end
+
+# Keep the original greet to avoid breaking any external references.
+greet() = print("Hello World!")
 
 end # module
*** End Patch
*** Begin Patch
*** Add File: test/runtests.jl
+using Test
+using Verlet
+
+@testset "ParticleSystem construction" begin
+    ps = ParticleSystem([0.0 0.0], [1.0 0.0], [1.0])
+    @test size(ps.positions) == (1,2)
+    @test size(ps.velocities) == (1,2)
+    @test ps.masses == [1.0]
+end
+
+@testset "Velocity Verlet free particle" begin
+    forces(r) = zeros(size(r))
+    ps = ParticleSystem([0.0 0.0], [1.0 0.0], [1.0])
+    velocity_verlet!(ps, forces, 0.1)
+    @test isapprox(ps.positions[1,1], 0.1; atol=1e-12)
+    @test isapprox(ps.velocities[1,1], 1.0; atol=1e-12)
+end
+
+@testset "Velocity Verlet harmonic oscillator" begin
+    forces(r) = -r
+    ps = ParticleSystem([1.0 0.0], [0.0 0.0], [1.0])
+    velocity_verlet!(ps, forces, 0.1)
+    @test ps.positions[1,1] < 1.0   # should move left
+end
*** End Patch

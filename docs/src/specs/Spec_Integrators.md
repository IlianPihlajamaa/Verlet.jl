# Spec: Module `Verlet.Integrators`

Purpose: Concrete `AbstractIntegrator` implementations that operate on `Verlet.Core.System` objects using the shared forcefield plumbing.

## Shared utilities
- `_ensure_forcefield(system::System)` validates that `system.forcefield` is attached; otherwise an `ArgumentError` is thrown.
- `_update_forces!(system::System)` delegates to `compute_all_forces!(system, system.forcefield)` and returns `system.forces`.
- Lightweight vector helpers (`_vecdot`, `_neg!`, `_copy!`) support conjugate-gradient updates.

## Velocity-Verlet
- `struct VelocityVerlet{T} <: AbstractIntegrator`
  - Fields: `dt::T`, `wrap::Bool` (wrap positions after each drift when `true`).
  - Constructor: `VelocityVerlet(dt::Real; wrap::Bool=false)`.
- `step!(integrator::VelocityVerlet, system::System)`
  - Calls `_update_forces!`, performs the usual half-kick → drift → optional wrap → force refresh → half-kick sequence, and updates velocities/positions in place.

## Conjugate-gradient minimiser
- `mutable struct ConjugateGradient{T,E} <: AbstractIntegrator`
  - Parameters: tolerance `tol`, initial step `alpha0`, minimum step `min_alpha`, Armijo slope coefficient `c1`, wrapping flag `wrap`, and user-supplied `energy::E` callable.
  - Holds a mutable `state::CGState` record with cached gradients, search direction, backup positions, latest energy, and stopping flags.
  - Constructor: `ConjugateGradient(energy; tol=1e-8, alpha0=1.0, min_alpha=1e-8, c1=1e-4, wrap=true)` promotes numeric types automatically.
- `step!` lazily initialises state buffers, computes forces/energy, performs a Polak–Ribière update with Armijo backtracking, and resets the search when slope conditions fail. On convergence or failed line search, the integrator sets `state.stop = true` so `stop_requested(integrator)` halts the outer `integrate!` loop.

## Langevin BAOAB
- `struct LangevinBAOAB{T,RNG} <: AbstractIntegrator`
  - Fields: `dt`, `gamma`, `temp`, `kB`, `wrap`, and `rng` (any `AbstractRNG`).
  - Constructor: `LangevinBAOAB(dt; γ, temp, kB=1.0, wrap=false, rng=Random.default_rng())` promotes numeric types.
- `step!` performs the BAOAB splitting: deterministic half-kick with current forces, half drift, Ornstein–Uhlenbeck velocity update using shared random normals, second half drift, and final half kick with refreshed forces. Optional wrapping applies after each drift.

## Constrained Langevin BAOAB
- `struct LangevinBAOABConstrained{T,RNG} <: AbstractIntegrator`
  - Extends the unconstrained variant with a `constraints::DistanceConstraints` field (from `Verlet.Constraints`).
  - Constructor mirrors the unconstrained signature: `LangevinBAOABConstrained(dt, constraints; γ, temp, kB=1.0, wrap=false, rng=Random.default_rng())`.
- `step!` interleaves SHAKE/RATTLE projections after each sub-step so both positions and velocities satisfy the supplied distance constraints.

## Behaviour
- All integrators rely on `_update_forces!`, ensuring neighbour lists are prepared through `Neighbors.prepare_neighbors!` before forces are consumed.
- `integrate!` is provided by `Verlet.Core`; callbacks returning `false` or `stop_requested(integrator)` terminating early are honoured for every integrator.
- Time-stepping mutates the supplied `System` in place; no copies of particle arrays are made unless an integrator explicitly caches them (e.g. conjugate gradient backups).

## Example
```julia
using Verlet, StaticArrays

box = CubicBox(5.0)
R = [SVector{3}(randn(), randn(), randn()) for _ in 1:8]; wrap_positions!(R, box)
sys = System(R, fill(SVector{3}(0.0, 0.0, 0.0), 8), fill(SVector{3}(0.0, 0.0, 0.0), 8), ones(8), box,
             ones(Int, 8), Dict(1 => :A); forcefield=ForceField(()))

vv = VelocityVerlet(5e-3)
integrate!(vv, sys, 10)
```

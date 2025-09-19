# Spec: Module `Verlet.Core`

Purpose: Core particle system types, math utilities, simple integrators, and force orchestration.

## Constants

- `T_Float = Float64`: Default floating type for numeric fields.
- `T_Int = Int64`: Default integer type for indices/types.
- `Dims = 3`: Default spatial dimensionality (APIs accept any `D`).

## Boxes

- `abstract type AbstractBox{T<:AbstractFloat}`: Base for boundary conditions.
- `struct CubicBox{T<:AbstractFloat} <: AbstractBox{T}`
  - Fields: `L::T` (box side length).
  - `box_length(box::CubicBox) -> T`
  - `minimum_image(Δ::AbstractVector, box::CubicBox)` and `minimum_image(Δ, L)`
    - Wrap displacement to (−L/2, L/2] componentwise.
  - `wrap_positions!(R::Vector{SVector{D,T}}, box::CubicBox{T})`
    - In-place wrap of positions to primary cell using minimum image.

## System and Particles

- `struct System{T<:AbstractFloat,IT<:Integer,Dims}`
  - Fields:
    - `positions::Vector{SVector{Dims,T}}`
    - `velocities::Vector{SVector{Dims,T}}`
    - `forces::Vector{SVector{Dims,T}}`
    - `masses::Vector{T}`
    - `box::AbstractBox{T}`
    - `types::Vector{IT}`
    - `type_names::Dict{IT,Symbol}`
    - `natoms::IT` (derived at construction)
    - `specific_potentials::Tuple` (bonded interactions)
    - `forcefield`: object supplying `compute_forces!` (typically a `ForceField`) or `nothing`
  - Invariants:
    - All particle arrays have equal length `natoms`.
    - `positions[i]`, `velocities[i]`, `forces[i]` are length-`Dims` `SVector`s.
- `natoms(sys::System) -> Integer`
- `natomtypes(sys::System) -> Int`
- `kinetic_energy(sys::System{T}) -> T`

## Abstract Potential Tags

- `AbstractPotentialPair`, `AbstractPairPotential`, `AbstractBondPotential`, `AbstractAnglePotential`, `AbstractDihedralPotential`, `AbstractImproperPotential`.
  - Marker types used by potentials for dispatch; no fields or behavior here.

## Neighbor Types (tags)

- `abstract type AbstractNeighborList end`
  - Concrete neighbor structures live in `Neighbors` and are used by Core APIs via dispatch.

## Forces and ForceField

- `struct ForceField{ForcesTuple}`
  - `layers::ForcesTuple` where each layer is a concrete potential (e.g., `LennardJones`, `Coulomb`).
- Generic multimethod hooks:
  - `compute_forces!(pot, sys::System)` must be implemented by each potential.
- `compute_all_forces!(sys::System, ff::ForceField)`
    - Resets `sys.forces` to zero element and accumulates contributions from each layer in `ff.layers`, followed by bonded `sys.specific_potentials`.
    - Requires per-layer neighbor lists to be prepared beforehand (e.g., via `Neighbors.build_all_neighbors!`).
  - `compute_all_forces!(sys::System)`
    - Uses the `sys.forcefield`; throws if none is attached.
- `compute_potential_energy(sys::System, ff::ForceField)` / `compute_potential_energy(sys::System)`
  - Ensures neighbor information is up to date and returns the scalar potential energy from the `ForceField` layers and any bonded interactions.

- `AbstractIntegrator`
  - Marker supertype for all time integration / minimisation schemes.
- `integrate!(integrator::AbstractIntegrator, system::System, nsteps; callback=nothing, kwargs...)`
  - Drives the integration loop: rebuilds neighbor lists via `rebuild_neighbors!`, calls `step!(integrator, system)`, and invokes the optional callback. Callbacks receive `(system, step, integrator)` and may return `false` to halt early. The attached `ForceField` stores and maintains its own master neighbour list when pair potentials are present.
- `step!(integrator::AbstractIntegrator, system::System)`
  - Each concrete integrator implements one-step updates; return value is ignored.
- `rebuild_neighbors!(system, args...; kwargs...)`
  - Hook for neighbor maintenance. `Verlet.Neighbors` provides methods for `MasterNeighborList`.
- `stop_requested(integrator::AbstractIntegrator)`
  - Integrators may signal convergence/termination via this query (default `false`).
- `potential_energy(system::System{T}, forces::Function) -> T`
  - Expects `forces(R; return_potential=true) => (F, U)`; errors if unsupported.

## Performance & Complexity

- All per-particle loops are `O(N)`. Force complexity depends on neighbor list size and potential implementations.
- `wrap_positions!` and `minimum_image` are allocation-free for `SVector` inputs.

## Error Handling

- `System` constructor asserts matching array sizes.
- `potential_energy` throws if `forces` does not accept the `return_potential=true` keyword.

## Example

```julia
using Verlet.Core, Verlet.Integrators, StaticArrays
box = CubicBox(10.0)
R = [@SVector randn(3) for _ in 1:4]; wrap_positions!(R, box)
sys = System(R, fill(@SVector zeros(3), 4), fill(@SVector zeros(3), 4), ones(4), box, ones(Int,4), Dict(1=>:A);
           forcefield=ForceField(()))
vv = VelocityVerlet(0.001)
integrate!(vv, sys, 1)
```

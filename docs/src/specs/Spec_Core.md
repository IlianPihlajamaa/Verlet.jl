# Spec: Module `Verlet.Core`

Purpose: Foundational particle-system types, periodic-box utilities, abstract potential tags, and shared integration hooks.

## Constants
- `T_Float = Float64`: default floating-point element type used by convenience constructors.
- `T_Int = Int64`: default integer type for particle indices and type labels.
- `Dims = 3`: default spatial dimensionality; most APIs work for any `D` through dispatch.

## Boxes
- `abstract type AbstractBox{T<:AbstractFloat}`: supertype for boundary-condition containers.
- `struct CubicBox{T<:AbstractFloat} <: AbstractBox{T}`
  - Field: `L::T` (edge length of the cubic cell).
  - `box_length(box::CubicBox) -> T` returns the stored length.
  - `minimum_image(Δ::AbstractVector, box::CubicBox)` plus specialised `SVector` methods keep displacements in (−L/2, L/2] componentwise without allocations.
  - `wrap_positions!(R::Vector{SVector{D,T}}, box::CubicBox{T})` wraps particle positions to the primary cell in place.

## System container
- `struct System{T<:AbstractFloat,IT<:Integer,Dims,BOX<:AbstractBox{T},FF}`
  - Fields: `positions`, `velocities`, `forces` :: `Vector{SVector{Dims,T}}`; `masses::Vector{T}`; `box::BOX`; `types::Vector{IT}`; `type_names::Dict{IT,Symbol}`; `natoms::IT`; `specific_potentials::Tuple`; `forcefield::FF` (typically a `Neighbors.ForceField` or `nothing`).
  - Constructor checks that all particle arrays share a common length `natoms`.
- `natoms(sys::System)` and `natomtypes(sys::System)` expose counts derived from stored fields.
- `kinetic_energy(sys::System{T}) -> T` computes `∑ₖ ½ mₖ ‖vₖ‖²` allocation-free.

## Potential & neighbour tags
- `AbstractPotential`, `AbstractPairPotential`, `AbstractBondPotential`, `AbstractAnglePotential`, `AbstractDihedralPotential`, `AbstractImproperPotential`, and `AbstractPotentialPair` provide dispatch markers for concrete force implementations.
- `abstract type AbstractNeighborList end` is a placeholder supertype; concrete lists live in `Verlet.Neighbors`.

## Force orchestration
- `compute_forces!(pot::AbstractPotential, sys::System)` must be specialised by each potential; the Core fallback throws.
- `prepare_neighbors!(ff, sys::System)` is a hook for neighbour-aware forcefields. The Core method returns `ff` unchanged; `Verlet.Neighbors` specialises it for its `ForceField` type.
- `compute_all_forces!(sys::System, ff)`
  - Fills `sys.forces` with zeros, calls `prepare_neighbors!`, then invokes `compute_forces!` across `ff.layers` followed by every entry in `sys.specific_potentials`.
  - Returns the mutated `System` to support chaining.
- `compute_all_forces!(sys::System)` delegates through `sys.forcefield`, raising an `ArgumentError` if none is attached.
- `compute_potential_energy(sys::System, ff)` mirrors the force routine but accumulates scalar energy contributions via `compute_potential_energy(pot, sys)` (Core’s default for unknown potentials is zero).
- `compute_potential_energy(sys::System)` again delegates through `sys.forcefield` with the same guard.
- `potential_energy(system::System{T}, forces::Function) -> T` expects `forces(system.positions; return_potential=true) => (F, U)` and throws a descriptive error when the callable lacks that protocol.

## Observables Interface

- `abstract type Observable end`: protocol marker defined in `Verlet.Core`.
- `observe!(obs::Observable, system::System, step::Integer)`
  - Implementations compute the value for `step` and return it (or a reusable scratch buffer).
  - Observables never own the persistent storage for their measurements; loggers manage retention and scheduling.
- `observed_quantity(obs)::Symbol` -> optional helper returning canonical name.
- `out_type(obs)::DataType` returns the type of an observation.
- Observables are registered with loggers in `Verlet.Loggers`; Core provides only the minimal interface so observables can be authored without depending on logger implementations.

## Integration interface
- `abstract type AbstractIntegrator end`: parent type for all integrators.
- `integrate!(integrator::AbstractIntegrator, system::System, nsteps::Integer; callback=nothing, neighbor_kwargs...)`
  - Validates `nsteps ≥ 0`; `nsteps == 0` returns immediately.
  - Calls `step!(integrator, system)` for each iteration and invokes the optional callback with `(system, step, integrator)`; if the callback returns `false` or `stop_requested(integrator)` is `true`, the loop exits early.
- `step!(integrator::AbstractIntegrator, system::System)` must be implemented by concrete integrators; the Core definition raises a `MethodError`.
- `rebuild_neighbors!` is declared with no methods in Core so that other modules (e.g. `Verlet.Neighbors`) can extend it.
- `maybe_rebuild(system::System, args...; kwargs...)` calls `rebuild_neighbors!` and returns the last positional argument, offering a lightweight helper for neighbour maintenance.
- `stop_requested(::AbstractIntegrator)` defaults to `false` but can be overridden by stateful integrators such as conjugate-gradient minimisation.

## Error handling & complexity
- `System` construction asserts consistent vector lengths; violations surface immediately.
- Force/energy helpers fall back to informative errors when a potential lacks the required overload.
- Core loops are `O(N)`; neighbour-related complexity is delegated to specialised modules.

## Example
```julia
using Verlet, StaticArrays

box = CubicBox(10.0)
R = [SVector{3}(randn(), randn(), randn()) for _ in 1:4]; wrap_positions!(R, box)
zeros3() = SVector{3}(0.0, 0.0, 0.0)

sys = System(R, fill(zeros3(), 4), fill(zeros3(), 4), ones(4), box,
             ones(Int, 4), Dict(1 => :A); forcefield=ForceField(()))

vv = VelocityVerlet(1e-3)
integrate!(vv, sys, 1)
```

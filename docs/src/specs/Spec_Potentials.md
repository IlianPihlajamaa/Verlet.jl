# Spec: Module `Verlet.Potentials`

Purpose: Concrete interatomic potentials, pair-parameter tables, and bonded interaction definitions that plug into the Core/Neighbors force orchestration.

## Pair-parameter infrastructure
- `struct LJPair{T<:Number} <: AbstractPotentialPair`
  - Fields: `ε::T`, `σ::T`, `rc::T` (cutoff radius).
- `struct CoulPair{T<:Number} <: AbstractPotentialPair`
  - Fields: `q1q2::T` (charge product) and `rc::T`.
- `struct PairTable{F<:AbstractPotentialPair}` holds a dense `Matrix{F}` mapping `(type_i, type_j)` to parameter objects.

## Lennard-Jones potential
- `mutable struct LennardJones{IntT<:Integer, T<:AbstractPotentialPair, T_Float} <: AbstractPairPotential`
  - Fields: `params::PairTable{T}`, `exclusions::Vector{Tuple{IntT, IntT}}`, `neighborlist::PotentialNeighborList{T, IntT}`, `skin::T_Float`.
  - Constructor: `LennardJones(params::PairTable{T}, exclusions, skin)` initialises `neighborlist` with `PotentialNeighborList(eltype(params.table))`.
- `compute_forces!(pot::LennardJones, sys::System)`
  - Iterates `pot.neighborlist.neighbors`, applies minimum-image displacements, and accumulates the 12-6 force when `r^2 < rc^2`.
- `compute_potential_energy(pot::LennardJones, sys::System)` sums the standard LJ energy for neighbours within the cutoff.

## Coulomb potential
- `mutable struct Coulomb{IntT<:Integer, T<:AbstractPotentialPair, T_Float} <: AbstractPairPotential`
  - Fields: `params::PairTable{T}`, `exclusions::Vector{Tuple{IntT, IntT}}`, `neighborlist::PotentialNeighborList{T}`, `skin::T_Float` (matching the current struct definition).
  - Constructor: `Coulomb(params::PairTable{T}, exclusions, skin)` also uses `PotentialNeighborList(eltype(params.table))`.
- `compute_forces!(pot::Coulomb, sys::System)` computes `f = (q₁q₂ / r³) Δ` for neighbours within `rc`.
- `compute_potential_energy(pot::Coulomb, sys::System)` accumulates `q₁q₂ / r` for valid neighbours, skipping zero-distance pairs.

## Bonded interactions
- Parameter types:
  - `struct HarmonicBond{T} <: AbstractBondPotential` (`k`, `r0`).
  - `struct HarmonicAngle{T} <: AbstractAnglePotential` (`k`, `θ0`).
  - `struct PeriodicDihedral{T} <: AbstractDihedralPotential` (`k`, `n::Int`, `ϕ0`).
- Instances binding particle indices to parameters:
  - `struct Bond{T<:AbstractBondPotential}` with `i`, `j`, `potential::T`.
  - `struct Angle{T<:AbstractAnglePotential}` with `i`, `j`, `k`, `potential::T`.
  - `struct Dihedral{T<:AbstractDihedralPotential}` with `i`, `j`, `k`, `l`, `potential::T`.
- Force/energy implementations:
  - `compute_forces!(bond::Bond, system)` applies the harmonic bond force aligned with the bond axis; `compute_potential_energy` returns `½k(r − r₀)²`.
  - `compute_forces!(angle::Angle, system)` implements the OpenMM-style harmonic angle force; `compute_potential_energy` uses `½k(θ − θ₀)²`.
  - `compute_forces!(dihedral::Dihedral, system::System{T,IT,3})` evaluates the periodic torsion using 3D cross products; the 2D/ND fallback throws a `DomainError` indicating the dimensionality requirement. Matching `compute_potential_energy` returns `k(1 + cos(nϕ − ϕ₀))` in 3D and zero otherwise.

## Integration with `Verlet.Neighbors`
- Pair potentials expect to be part of a `Neighbors.ForceField`. During `compute_all_forces!` the specialised `prepare_neighbors!` builds or updates a shared `MasterNeighborList` and re-populates each `pot.neighborlist` according to `params` and `skin`.
- Bonded interactions should be provided through `system.specific_potentials` so they are included automatically after the pairwise layers.
- Exclusions are handled by `Neighbors.is_excluded`, which checks tuple membership in `pot.exclusions`.

## Behaviour & invariants
- Pair potential kernels are allocation-free over `SVector` positions.
- Neighbour lists are cleared and rebuilt on demand; capacity is reused between builds.
- `skin` expands the accepted neighbour radius to `(rc + skin)` when filtering candidates from the master list.

## Example
```julia
using Verlet, StaticArrays

box = CubicBox(6.0)
R = [SVector{3}(randn(), randn(), randn()) for _ in 1:8]; wrap_positions!(R, box)
sys = System(R, fill(SVector{3}(0.0, 0.0, 0.0), 8), fill(SVector{3}(0.0, 0.0, 0.0), 8), ones(8), box,
             ones(Int, 8), Dict(1 => :A))

params = PairTable(fill(LJPair(1.0, 1.0, 2.5), (1, 1)))
lj = LennardJones(params, Tuple{Int,Int}[], 0.4)
ff = ForceField((lj,))
compute_all_forces!(sys, ff)
```

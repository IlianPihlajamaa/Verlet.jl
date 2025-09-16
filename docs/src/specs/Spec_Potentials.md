# Spec: Module `Verlet.Potentials`

Purpose: Define concrete interatomic potentials and per-pair parameter tables; integrate with neighbor lists for efficient force evaluation.

## Pair parameter structs

- `struct LJPair{T<:Number} <: AbstractPotentialPair`
  - `ε::T`, `σ::T`, `rc::T` (cutoff).
- `struct CoulPair{T<:Number} <: AbstractPotentialPair`
  - `q1q2::T` (charge product), `rc::T`.
- `struct PairTable{F<:AbstractPotentialPair}`
  - `table::Matrix{F}` mapping `(type_i, type_j)` to parameter objects.

## Concrete pair potentials

- `mutable struct LennardJones{IntT<:Integer, T<:AbstractPotentialPair, T_Float} <: AbstractPairPotential`
  - Fields: `params::PairTable{T}`, `exclusions::Vector{Tuple{IntT,IntT}}`, `neighbors::PotentialNeighborList{T}`, `skin::T_Float`.
  - Constructor: `LennardJones(params::PairTable{T}, exclusions, skin)` initializes empty neighbor list.
  - `compute_forces!(pot::LennardJones, sys::System)`
    - For each neighbor `(i,j,p)`, computes minimum-image displacement `Δ`, `r2 = Δ⋅Δ`.
    - If `r2 < p.rc^2`, apply LJ 12-6 force: `f = 24ε(2(σ^12/r^13) - (σ^6/r^7)) Δ/‖Δ‖` (implemented in r2 form for efficiency) and accumulate on `sys.forces[i]`/`[j]`.
- `mutable struct Coulomb{IntT<:Integer, T<:AbstractPotentialPair, T_Float} <: AbstractPairPotential`
  - Fields analogous to `LennardJones` but with `CoulPair` parameters.
  - `compute_forces!(pot::Coulomb, sys::System)`
    - For each neighbor within cutoff, `f = (q1q2 / r^3) Δ` (with minimum-image) and accumulate.

## Bonded potentials and instances

- Potential parameters
  - `struct HarmonicBond{T} <: AbstractBondPotential` with `k`, `r0`.
  - `struct HarmonicAngle{T} <: AbstractAnglePotential` with `k`, `θ0`.
  - `struct PeriodicDihedral{T} <: AbstractDihedralPotential` with `k`, `n::Int`, `ϕ0`.
- Interaction instances (bind particles to parameters)
  - `struct Bond{T<:AbstractBondPotential}` with particles `(i,j)` and `potential::T`.
  - `struct Angle{T<:AbstractAnglePotential}` with `(i,j,k)` and `potential::T`.
  - `struct Dihedral{T<:AbstractDihedralPotential}` with `(i,j,k,l)` and `potential::T`.
- Force implementations
  - `compute_forces!(bond::Bond, system)` harmonic stretch along bond axis.
  - `compute_forces!(angle::Angle, system)` harmonic deviation of angle at `j`.
  - `compute_forces!(dihedral::Dihedral, system)` periodic torsion per OpenMM formulas.

## Integration with `Core.ForceField`

- Pair potentials are included in `ForceField((lj, coul, ...))` and must have neighbors prebuilt via `Neighbors.build_all_neighbors!`.
- Bonded interactions should be placed in `system.specific_potentials` for automatic inclusion by `compute_all_forces!`.

## Exclusions

- Each pair potential carries `exclusions::Vector{Tuple{Int,Int}}` checked by `Neighbors.is_excluded` during neighbor construction.

## Invariants & Performance

- `neighbors` lists are cleared and rebuilt; structure capacity may be retained to avoid allocations.
- Force kernels run allocation-free and scale with neighbor count.

## Examples

```julia
using Verlet, StaticArrays
box = CubicBox(8.0)
R = [@SVector randn(3) for _ in 1:64]; wrap_positions!(R, box)
sys = System(R, fill(@SVector zeros(3),64), fill(@SVector zeros(3),64), ones(64), box, ones(Int,64), Dict(1=>:A))

params = PairTable(fill(LJPair(1.0, 1.0, 2.5), (1,1)))
lj = LennardJones(params, Tuple{Int,Int}[], 0.5)
ff = Verlet.Neighbors.ForceField((lj,))
master = Verlet.Neighbors.MasterNeighborList(sys; cutoff=2.5, skin=0.5)
Verlet.Neighbors.build_all_neighbors!(master, ff, sys)
Verlet.Neighbors.compute_all_forces!(sys, ff)
```

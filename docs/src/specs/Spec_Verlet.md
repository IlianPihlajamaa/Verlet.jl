# Spec: Module `Verlet`

Purpose: Assemble the submodules (`Core`, `Neighbors`, `Potentials`, `Constraints`, `Integrators`, `Thermostats`, and the WIP `Electrostatics`) into a cohesive user-facing API.

## Re-exports
- Core & integration: `System`, `natoms`, `natomtypes`, `AbstractBox`, `CubicBox`, `minimum_image`, `wrap_positions!`, `box_length`, `AbstractIntegrator`, `integrate!`, `VelocityVerlet`, `ConjugateGradient`, `LangevinBAOAB`, `LangevinBAOABConstrained`, `potential_energy`, `compute_potential_energy`, `kinetic_energy`, `maybe_rebuild`, `T_Float`, `T_Int`, `Dims`.
- Neighbour infrastructure: `NeighborPair`, `PotentialNeighborList`, `MasterNeighborList`, `ForceField`, `build_master_neighborlist!`, `build_cellgrid`, `rebin!`, `brute_force_pairs`.
- Potential markers & concretes: `AbstractPotentialPair`, `AbstractPairPotential`, `AbstractBondPotential`, `AbstractAnglePotential`, `AbstractDihedralPotential`, `AbstractImproperPotential`, `LennardJones`, `Coulomb`, `LJPair`, `CoulPair`, `PairTable`, `Bond`, `Angle`, `Dihedral`, `HarmonicBond`, `HarmonicAngle`, `PeriodicDihedral`.
- Constraints: `DistanceConstraints`, `apply_shake!`, `apply_rattle!`, `velocity_verlet_shake_rattle!`, `remove_com_motion!`, `constraint_residuals`.
- Thermostats: `degrees_of_freedom`, `instantaneous_temperature`, `velocity_rescale!`.

## Semantics
- `using Verlet` is sufficient for typical workflows; advanced users can still access submodules directly (e.g. `Verlet.Neighbors.ForceField`).
- Re-exported names follow the project’s semantic-versioning policy; non-exported internals are considered unstable.

## Example
```julia
using Verlet, StaticArrays

box = CubicBox(10.0)
R = [SVector{3}(randn(), randn(), randn()) for _ in 1:64]; wrap_positions!(R, box)
sys = System(R, fill(SVector{3}(0.0, 0.0, 0.0), 64), fill(SVector{3}(0.0, 0.0, 0.0), 64), ones(64), box,
             ones(Int, 64), Dict(1 => :A))
ke = kinetic_energy(sys)
```

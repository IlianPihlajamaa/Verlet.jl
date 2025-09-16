# Spec: Module `Verlet`

- Purpose: Provide a cohesive top-level API by loading internal submodules and re-exporting core types and functions for users.
- Submodules: `Core`, `Neighbors`, `Potentials`, `Constraints`, `Thermostats` (and `Electrostatics` WIP).

## Re-exports (user-facing)

- Core: `System`, `natoms`, `natomtypes`, `AbstractBox`, `CubicBox`, `minimum_image`, `wrap_positions!`, `box_length`, `velocity_verlet!`, `integrate!`, `potential_energy`, `kinetic_energy`, `T_Float`, `T_Int`, `Dims`.
- Neighbor types: `NeighborPair`, `PotentialNeighborList`, `MasterNeighborList`.
- Potential types: `AbstractPotentialPair`, `AbstractPairPotential`, `AbstractBondPotential`, `AbstractAnglePotential`, `AbstractDihedralPotential`, `AbstractImproperPotential`.
- Concrete potentials: `LennardJones`, `Coulomb`, `LJPair`, `CoulPair`, `PairTable`, `Bond`, `Angle`, `Dihedral`, `HarmonicBond`, `HarmonicAngle`, `PeriodicDihedral`.
- Neighbor building / forcefield: `build_master_neighborlist!`, `build_cellgrid`, `rebin!`, `ForceField`, `brute_force_pairs`.
- Constraints: `DistanceConstraints`, `apply_shake!`, `apply_rattle!`, `velocity_verlet_shake_rattle!`, `remove_com_motion!`, `constraint_residuals`.
- Thermostats: `degrees_of_freedom`, `instantaneous_temperature`, `velocity_rescale!`, `langevin_baoab!`, `langevin_baoab_constrained!`.

## Semantics

- Importing `Verlet` is sufficient for most workflows; the above names are available from the top-level module.
- Submodules can still be accessed directly (e.g., `Verlet.Neighbors.ForceField`).

## Versioning & Stability

- Public names re-exported here follow semantic versioning; removals or breaking signature changes increment the major version.
- Internals not re-exported are considered unstable.

## Example

```julia
using Verlet, StaticArrays
box = CubicBox(10.0)
R = [@SVector randn(3) for _ in 1:64]; wrap_positions!(R, box)
sys = System(R, fill(@SVector zeros(3), 64), fill(@SVector zeros(3), 64), ones(64), box, ones(Int,64), Dict(1=>:A))
KE = kinetic_energy(sys)
```

# Bonded Interactions

Bonded terms capture interactions tied to specific atom tuples—bonds, angles,
and dihedrals. They complement pair potentials and are evaluated automatically
once added to the `specific_potentials` tuple on a `System`.

## Available interaction types

- `HarmonicBond`, `HarmonicAngle`, and `PeriodicDihedral` define the parameter families for bonds, angles, and torsions.
- `Bond`, `Angle`, and `Dihedral` bind those parameters to explicit atom indices.

## Example: triatomic fragment

```@example bonded
using StaticArrays, Verlet

zero3 = SVector{3}(0.0, 0.0, 0.0)
positions = [zero3,
             SVector{3}(1.0, 0.0, 0.0),
             SVector{3}(1.0, 1.0, 0.0)]
velocities = fill(zero3, 3)
forces = fill(zero3, 3)
masses = ones(3)
box = CubicBox(10.0)
types = [1, 1, 1]
type_names = Dict(1 => :A)

bond = Bond(1, 2, HarmonicBond(200.0, 1.0))
angle = Angle(1, 2, 3, HarmonicAngle(50.0, deg2rad(109.5)))

ff = ForceField(())
data = System(positions, velocities, forces, masses, box, types, type_names;
              forcefield = ff, specific_potentials = (bond, angle))

Verlet.Core.compute_all_forces!(data)
(data.forces[1], data.forces[2])
```

## Workflow tips

- Supply bonded terms via the `specific_potentials` keyword when constructing a
  `System`, or mutate the tuple after construction (`sys = sys`.with etc.).
- Combine as many interactions as required: the tuple is iterated in order, so
  you can include custom objects that also implement `compute_forces!`.
- Use bonded terms alongside neighbour-managed pair potentials by attaching a
  `ForceField` and setting `specific_potentials` simultaneously.

## Debugging

- If a bonded interaction fails to converge or produces unexpected forces,
  inspect the angle/dihedral geometry directly—ill-conditioned setups (nearly
  colinear atoms) can amplify numerical error.
- `compute_potential_energy(bond, system)` and friends return the contribution
  from an individual interaction; sum them to decompose total energy.

For constrained alternatives (fixed bond lengths), see
[Tutorial 3 · Constraints in Practice](../tutorials/constraints.md).

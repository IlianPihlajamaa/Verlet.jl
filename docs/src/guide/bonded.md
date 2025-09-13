# Bonded Potentials

Verlet.jl supports a variety of bonded interactions, which are forces that apply to specific groups of atoms (e.g., bonds, angles, dihedrals). These are specified by passing a tuple of interaction objects to the `System` constructor via the `specific_potentials` keyword argument.

## Available Potentials

### Harmonic Bond

The harmonic bond potential is given by:

```math
U(r) = \frac{1}{2} k (r - r_0)^2
```

It is defined in `Verlet.jl` as `HarmonicBond(k, r0)`. To apply this potential to a pair of atoms, you create a `Bond` object:

```julia
bond_potential = HarmonicBond(1000.0, 0.15) # k in kJ/mol/nm^2, r0 in nm
bond = Bond(1, 2, bond_potential) # Applies to atoms 1 and 2
```

### Harmonic Angle

The harmonic angle potential is given by:

```math
U(\theta) = \frac{1}{2} k (\theta - \theta_0)^2
```

It is defined as `HarmonicAngle(k, θ0)`. To apply it to a triplet of atoms, you create an `Angle` object:

```julia
angle_potential = HarmonicAngle(100.0, deg2rad(109.5)) # k in kJ/mol/rad^2, θ0 in radians
angle = Angle(1, 2, 3, angle_potential) # Applies to atoms 1-2-3 (2 is central)
```

### Periodic Dihedral

The periodic dihedral (or torsion) potential is given by:

```math
U(\phi) = k (1 + \cos(n\phi - \phi_0))
```

It is defined as `PeriodicDihedral(k, n, ϕ0)`. To apply it to a quadruplet of atoms, you create a `Dihedral` object:

```julia
dihedral_potential = PeriodicDihedral(10.0, 3, 0.0) # k in kJ/mol, n is integer, ϕ0 in radians
dihedral = Dihedral(1, 2, 3, 4, dihedral_potential)
```

## Example Usage

To use these in a simulation, group them into a tuple and pass them to the `System`:

```julia
bond1 = Bond(1, 2, HarmonicBond(1000.0, 0.15))
angle1 = Angle(1, 2, 3, HarmonicAngle(100.0, deg2rad(109.5)))

specific_interactions = (bond1, angle1)

sys = System(
    ...,
    specific_potentials=specific_interactions
)
```

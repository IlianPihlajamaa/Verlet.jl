# API Reference

```@docs
Verlet
ParticleSystem
CubicBox
velocity_verlet!
kinetic_energy
potential_energy
minimum_image!
lj_forces
```

## Neighbor Lists

```@docs
NeighborList
build_neighborlist
maybe_rebuild!
max_displacement_since_build
wrap_positions!
```

## Cell Grid + O(N) Builder (Half Lists)

```@docs
CellGrid
build_cellgrid
rebin!
build_neighborlist_cells
```

## Thermostats

```@docs
langevin_baoab!
instantaneous_temperature
degrees_of_freedom
velocity_rescale!
```

## Notes

Some internal helper functions (`_box_length`, `_linear_index`, etc.) are not exported
and intentionally excluded from this API reference.

## Constraints

```@docs
DistanceConstraints
velocity_verlet_shake_rattle!
apply_shake!
apply_rattle!
remove_com_motion!
```

## Numerical Notes on Constraints

Constraint solvers may fail to converge if `tol` is set too small or if bond networks are
ill-conditioned. Use `maxiter` to cap the iterations, and consider loosening `tol` if you see convergence errors.

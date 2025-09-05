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

## Notes

Some internal helper functions (`_box_length`, `_linear_index`, etc.) are not exported
and intentionally excluded from this API reference.

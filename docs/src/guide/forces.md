# Forces & Potentials

A **force function** takes positions `r::Matrix{Float64}` of shape `(N×D)` and returns a force matrix of the same shape.

```@example forces
# Free particle (no forces)
forces_free(r) = zeros(size(r))

# Linear spring to the origin (Hooke's law, k = 1)
forces_ho(r) = -r

# With potential-energy support:
function forces_ho_with_U(r; return_potential=false)
    F = -r
    return return_potential ? (F, 0.5 * sum(abs2, r)) : F
end
nothing
# Forces & Potentials
```
This page collects simple force examples and conventions used by `Verlet.jl`.

## Lennard–Jones with/without Neighbor Lists

```@example forces
using Verlet, LinearAlgebra
box = CubicBox(10.0)
R = randn(64,3); wrap_positions!(R, box) 
# Brute force (O(N²))
F_bf, U_bf = lj_forces(R, box; rcut=2.5, return_potential=true) 
# Classic symmetric neighbor list (O(N) per step, O(N²) build)
nlsym = build_neighborlist(R, box; cutoff=2.5, skin=0.4)
F_sym, U_sym = lj_forces(R, box, nlsym; rcut=2.5, return_potential=true) 
# Cell-based half neighbor list (O(N) build, O(N) per step)
grid = build_cellgrid(R, box; cell_size=2.9)
nlhalf = build_neighborlist_cells(R, box; cutoff=2.5, skin=0.4, grid=grid)
F_half, U_half = lj_forces(R, box, nlhalf; rcut=2.5, return_potential=true) 
(norm(F_bf - F_sym), norm(F_bf - F_half))
```

## Custom forces with potential

To make `potential_energy` work, return `(F, U)` when the keyword
`return_potential=true` is provided:

```@example forces2
using Verlet
function ho_forces(r; return_potential=false)
    F = -r      # k = 1
    return return_potential ? (F, 0.5 * sum(abs2, r)) : F
end
ps = ParticleSystem([1.0 0.0], [0.0 0.0], [1.0])
E = potential_energy(ps, ho_forces)
E
```

## Performance notes

- Prefer **half lists** for LJ when possible; they reduce memory and branches.
- With **O(N) builds**, you can lower `skin` (e.g., 0.2–0.3) for tighter forces.
- Keep arrays as `Float64` to minimize drift.
# Forces & Potentials

```@example forces
using StaticArrays
# Free particle (no forces)
forces_free(R) = [@SVector zeros(length(R[1])) for _ in R]

# Linear spring to the origin (Hooke's law, k = 1)
forces_ho(R) = [ -r for r in R ]

# With potential-energy support:
function forces_ho_with_U(R; return_potential=false)
    F = [ -r for r in R ]
    U = 0.5 * sum(norm(r)^2 for r in R)
    return return_potential ? (F, U) : F
end
nothing
# Forces & Potentials
```
## Lennard–Jones with/without Neighbor Lists

```@example forces
using Verlet, StaticArrays, LinearAlgebra

box = CubicBox(10.0)
R = [@SVector(randn(3)) for _ in 1:64]; wrap_positions!(R, box)
# Brute force (O(N²))
F_bf, U_bf = lj_forces(R, box; rcut=2.5, return_potential=true)
# Classic symmetric neighbor list (O(N) per step, O(N²) build)
nlsym = build_neighborlist(R, box; cutoff=2.5, skin=0.4)
F_sym, U_sym = lj_forces(R, box, nlsym; rcut=2.5, return_potential=true)
# Cell-based half neighbor list (O(N) build, O(N) per step)
grid = build_cellgrid(R, box; cell_size=2.9)
nlhalf = build_neighborlist_cells(R, box; cutoff=2.5, skin=0.4, grid=grid)
F_half, U_half = lj_forces(R, box, nlhalf; rcut=2.5, return_potential=true)
err_sym = sum(norm(F_bf[i] - F_sym[i]) for i in eachindex(F_sym))
err_half = sum(norm(F_bf[i] - F_half[i]) for i in eachindex(F_half)) 
(err_sym, err_half)
```

## Custom forces with potential

To make `potential_energy` work, return `(F, U)` when the keyword
`return_potential=true` is provided:

```@example forces2
using Verlet, StaticArrays, LinearAlgebra
function ho_forces(R; return_potential=false)
    F = [ -r for r in R ]      # k = 1
    U = 0.5 * sum(norm(r)^2 for r in R)
    return return_potential ? (F, U) : F
end
positions = [@SVector [1.0, 0.0]]
velocities = [@SVector [0.0, 0.0]]
masses = [1.0]
ps = ParticleSystem(positions, velocities, masses)
E = potential_energy(ps, ho_forces)
E
```

## Performance notes

- Prefer **half lists** for LJ when possible; they reduce memory and branches.
- With **O(N) builds**, you can lower `skin` (e.g., 0.2–0.3) for tighter forces.
- Use `Vector{SVector}` for all positions, velocities, and forces for best performance.
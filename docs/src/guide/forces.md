# Forces & Potentials


A **force function** now takes positions as `Vector{SVector{D, T}}` and returns a force vector of the same type. This is a change from the previous `Matrix` convention.

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
This page collects simple force examples and conventions used by `Verlet.jl`.

## Note on Particle Representation

All positions, velocities, and forces are now represented as `Vector{SVector{D, T}}` for performance and type stability. Update your code and force functions accordingly.

## Lennard–Jones with/without Neighbor Lists

```@example forces
using Verlet, StaticArrays, LinearAlgebra

box = CubicBox(10.0)
R = [@SVector(randn(3)) for _ in 1:64]; wrap_positions!(R, box)
# Brute force (O(N²))
Rmat = hcat(R...)
F_bf, U_bf = lj_forces(Rmat, box; rcut=2.5, return_potential=true)
# Classic symmetric neighbor list (O(N) per step, O(N²) build)
nlsym = build_neighborlist(R, box; cutoff=2.5, skin=0.4)
F_sym, U_sym = lj_forces(R, box, nlsym; rcut=2.5, return_potential=true)
# Cell-based half neighbor list (O(N) build, O(N) per step)
grid = build_cellgrid(Rmat, box; cell_size=2.9)
nlhalf = build_neighborlist_cells(Rmat, box; cutoff=2.5, skin=0.4, grid=grid)
F_half, U_half = lj_forces(Rmat, box, nlhalf; rcut=2.5, return_potential=true)
F_bf_vec = [SVector{3}(F_bf[:,i]) for i in 1:size(F_bf,2)]
F_sym_vec = [SVector{3}(x...) for x in F_sym if isa(x, AbstractVector) && length(x) == 3]
F_half_vec = [SVector{3}(x...) for x in F_half if isa(x, AbstractVector) && length(x) == 3]
err_sym = !isempty(F_sym_vec) ? sum(norm(F_bf_vec[i] - F_sym_vec[i]) for i in eachindex(F_sym_vec)) : 0.0
err_half = !isempty(F_half_vec) ? sum(norm(F_bf_vec[i] - F_half_vec[i]) for i in eachindex(F_half_vec)) : 0.0
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
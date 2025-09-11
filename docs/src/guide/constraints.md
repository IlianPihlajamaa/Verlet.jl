# Constrained Dynamics (SHAKE/RATTLE)

`@id constraints-guide
`

Molecular simulations often require certain bond lengths (e.g. X–H bonds in water) to remain fixed.
This enables **larger stable timesteps** and enforces realistic rigid-body structures.

Verlet.jl provides support for **holonomic distance constraints** using the
classical **SHAKE** (positions) and **RATTLE** (velocities) algorithms.

## Defining Constraints

Use [`DistanceConstraints`](@ref) to define a set of pairwise distance constraints:


```@example constraints
using Verlet, StaticArrays
# A diatomic molecule with target bond length 1.0
pairs   = [(1,2)]
lengths = [1.0]
cons = DistanceConstraints(pairs, lengths; tol=1e-10, maxiter=100)
```

Arguments:

- `pairs`: vector of `(i,j)` atom index pairs (1-based)
- `lengths`: vector of target distances
- `tol`: maximum squared violation tolerated (`|C_l|` units length²)
- `maxiter`: maximum SHAKE/RATTLE iterations per step
- `use_minimum_image`: apply minimum image convention under periodic boundaries

!!! tip
    For constraints across periodic boundaries, keep molecules whole and use
    `use_minimum_image=true`.

## Constrained Integrator

The [`velocity_verlet_shake_rattle!`](@ref) driver advances the system with constraints enforced:


```@example constraints
using Verlet, StaticArrays, LinearAlgebra
N, D = 2, 3
positions = [@SVector zeros(D) for _ in 1:N]
positions[2] = @SVector [1.2, 0.0, 0.0]   # initial bond slightly off
velocities = [@SVector zeros(D) for _ in 1:N]
forces = [@SVector zeros(D) for _ in 1:N]
masses = ones(N)
box = CubicBox(10.0)
types = ones(Int, N)
type_names = Dict(1 => :A)
sys = System(positions, velocities, forces, masses, box, types, type_names)
cons = DistanceConstraints([(1,2)], [1.0])
forces_func(R) = [@SVector zeros(D) for _ in R]  # no external forces
for step in 1:100
  velocity_verlet_shake_rattle!(sys, forces_func, 0.01, cons)
end
d = sys.positions[1] - sys.positions[2]
@show norm(d)  # ~1.0
```

This integrator:
1. Updates velocities (half step) and positions (drift).
2. Applies **SHAKE** projection to enforce bond lengths.
3. Recomputes forces.
4. Completes velocity update.
5. Applies **RATTLE** projection to enforce velocity constraints.

## Degrees of Freedom

Constraints reduce the effective number of degrees of freedom (DoF).
The [`degrees_of_freedom`](@ref) function accounts for:

- number of atoms × dimensions
- minus one per constraint
- minus dimensions if COM motion is removed


```@example constraints
using StaticArrays
N, D = 3, 3
positions = [@SVector zeros(D) for _ in 1:N]
velocities = [@SVector zeros(D) for _ in 1:N]
forces = [@SVector zeros(D) for _ in 1:N]
masses = ones(N)
box = CubicBox(10.0)
types = ones(Int, N)
type_names = Dict(1 => :A)
sys = System(positions, velocities, forces, masses, box, types, type_names)
cons = DistanceConstraints([(1,2)], [1.0])
dof = degrees_of_freedom(sys; constraints=cons, remove_com=true)
@show dof
```

Correct DoF is essential for unbiased temperature and pressure estimators.

## Removing Center-of-Mass Motion

Use [`remove_com_motion!`](@ref) to zero the mass-weighted center-of-mass velocity or position.
This prevents unphysical drift of the entire system.


```@example constraints
using StaticArrays
N, D = 3, 3
positions = [@SVector zeros(D) for _ in 1:N]
velocities = [@SVector ones(D) for _ in 1:N]
forces = [@SVector zeros(D) for _ in 1:N]
masses = [1.0, 2.0, 3.0]
box = CubicBox(10.0)
types = ones(Int, N)
type_names = Dict(1 => :A)
sys = System(positions, velocities, forces, masses, box, types, type_names)
remove_com_motion!(sys; which=:velocity)
# After removal, COM velocity should be ~0
Vcom = sum(sys.masses .* map(v -> v[1], sys.velocities)) / sum(sys.masses)
@show Vcom
```

## Performance Notes

- SHAKE/RATTLE typically converge in a few iterations for tree-like molecules.
- For rings or stiff networks, increase `maxiter` or relax `tol`.
- Always monitor constraint residuals if using larger timesteps.
- Thermostat steps that randomize velocities should be followed by `apply_rattle!`.

## Common Pitfalls

- Constraints assume well-defined molecular topology. If your system uses PBC,

  * ensure constrained atoms belong to the same molecule and do not cross cell
    boundaries unexpectedly.
  * A too-tight tolerance can lead to slow or failed convergence.
  * DoF reduction is essential: forgetting to pass `constraints` or `remove_com`
    to [`degrees_of_freedom`](@ref) will bias temperature estimates.

## Further Reading

- Ryckaert, Ciccotti & Berendsen (1977), *Numerical integration of the Cartesian
  equations of motion of a system with constraints: molecular dynamics of n-alkanes*,
  J. Comp. Phys. 23(3).
- Andersen (1983), *RATTLE: A "velocity" version of the SHAKE algorithm for
  molecular dynamics calculations*, J. Comp. Phys. 52(1).

These classical references describe the original SHAKE and RATTLE algorithms
implemented here.


## Note on Particle Representation

All positions, velocities, and forces are now represented as `Vector{SVector{D, T}}` for performance and type stability. Update your code and constraint definitions accordingly.

## See Also

- [`System`](@ref): container for positions, velocities, and masses.
- [`velocity_verlet!`](@ref): unconstrained Velocity-Verlet integrator.
- [`degrees_of_freedom`](@ref): count effective translational degrees of freedom.
- [`remove_com_motion!`](@ref): eliminate center-of-mass drift.

For thermostatting with constraints, project velocities with
[`apply_rattle!`](@ref) after randomization steps to remain on the constraint manifold.

## Next Steps

You can now combine constrained dynamics with neighbor lists, Lennard-Jones forces,
and thermostats. See the [Forces & Potentials](@ref) guide for force field setup.

---

**Summary:** SHAKE/RATTLE constraints in Verlet.jl let you simulate rigid bonds,
stabilize molecules, and safely increase integration timesteps. Use them with care,
monitor convergence, and adjust tolerances as needed.

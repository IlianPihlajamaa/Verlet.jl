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
```

To make [`potential_energy`](@ref) work, support the keyword protocol:

```@example forces
# Called internally as: forces(r; return_potential=true) => (F, U)
nothing
```

If your force does **not** provide `return_potential=true`, `potential_energy` will throw with a clear message.

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
````

To make [`potential_energy`](@ref) work, support the keyword protocol:

```@example forces
# Called internally as: forces(r; return_potential=true) => (F, U)
nothing
```

If your force does **not** provide `return_potential=true`, `potential_energy` will throw with a clear message.

---

# Periodic Boxes & Lennard–Jones Forces

This page documents the new periodic boundary support and the Lennard–Jones (LJ) force provider.

## Quick Start

```@example lj_quickstart
using Verlet

# Two particles in 3D under LJ with periodic boundaries
R = [0.9 0.0 0.0;
     0.0 0.0 0.0]
V = zeros(2, 3)
M = ones(2)
ps = ParticleSystem(copy(R), copy(V), M)

box = CubicBox(20.0)

# Wrap an LJ force provider with the package's potential-energy convention
forces(r; return_potential=false) = return_potential ?
    lj_forces(r, box; return_potential=true) :
    lj_forces(r, box; return_potential=false)

E0 = kinetic_energy(ps) + potential_energy(ps, forces)
velocity_verlet!(ps, forces, 1e-3)
E1 = kinetic_energy(ps) + potential_energy(ps, forces)

(E0, E1, ps.positions)
```

## API

### `CubicBox`

```julia
CubicBox(L::Real)
```

A minimal cubic periodic box with side length `L`.

### `minimum_image!`

```julia
minimum_image!(Δ::AbstractVector, box::CubicBox)
```

Applies the minimum-image convention in-place so each component of `Δ` lies in `(-L/2, L/2]`.

### `lj_forces`

```julia
lj_forces(positions::AbstractMatrix, box::CubicBox;
          ϵ::Real=1.0, σ::Real=1.0, rcut::Real=Inf,
          shift::Bool=false, return_potential::Bool=false)
```

Computes Lennard–Jones pair forces with minimum-image convention. Returns an `(N×D)` force matrix or `(F, U)` if `return_potential=true`.

* **Potential:** `U(r) = 4ϵ[(σ/r)^12 - (σ/r)^6]`
* **Force magnitude:** `f(r) = 24ϵ(2(σ/r)^12 - (σ/r)^6)/r`

Cutoff behavior:

* Pairs with `r ≤ rcut` are included.
* If `shift=true` and `rcut < Inf`, a constant `-U(rcut)` is applied so `U(rcut)=0`.

## Worked Examples

### Energy conservation

```@example lj_energy
using Verlet

N = 8; D = 3
R = randn(N, D) .+ 0.5
V = zeros(N, D)
M = ones(N)
ps = ParticleSystem(copy(R), copy(V), M)
box = CubicBox(25.0)

forces(r; return_potential=false) = return_potential ?
    lj_forces(r, box; σ=1.0, ϵ=1.0, rcut=2.5, shift=true, return_potential=true) :
    lj_forces(r, box; σ=1.0, ϵ=1.0, rcut=2.5, shift=true, return_potential=false)

E0 = kinetic_energy(ps) + potential_energy(ps, forces)
for _ in 1:100
    velocity_verlet!(ps, forces, 5e-3)
end
E1 = kinetic_energy(ps) + potential_energy(ps, forces)

(E0, E1)
```

### Minimum-image demonstration

```@example minimum_image
using Verlet

box = CubicBox(10.0)
Δ = [ 6.0, -6.0, 0.1 ]
minimum_image!(Δ, box)
Δ
```

## Performance Tips

* Prefer `Float64` positions/velocities to reduce drift.
* Use `rcut` with `shift=true` for faster simulations and well-behaved energies.
* This implementation is `O(N²)`; for `N ≳ 10³`, neighbor lists are recommended.

## Pitfalls

* **Singularity at r→0:** avoid overlapping particles.
* **Cutoff discontinuity:** with `shift=false` both potential and force jump; with `shift=true` only the force jumps.
* **Box length:** keep system roughly inside the box for clarity.

EOF

```
```

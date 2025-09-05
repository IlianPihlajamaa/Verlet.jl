# Numerics & Pitfalls  {#numerics}

## Timestep stability

- Velocity Verlet is symplectic and stable for sufficiently small `dt`.
- If you see exploding energies or trajectories, reduce `dt`.
- Start with a small `dt` (e.g. `1e-3` in your time units) and increase cautiously.

## Precision

Use `Float64` across positions, velocities, and masses to limit drift. Distances, forces, and energies are accumulated in `Float64` by default.

## Shapes & broadcasting

- Arrays are `(N×D)` with rows as particles.
- `masses` is length `N`. Broadcasting `forces(r) ./ masses` divides each row by its particle mass.

## Boundary conditions (not included)

No periodic or reflective boundaries are included. Implement them in your force or position update logic if needed.

## Potential energy protocol

- Provide `forces(r; return_potential=true) => (F, U)` to enable `potential_energy`.
- If unknown, prefer returning only `F` and skip potential-energy reporting.

## Units

All quantities are unitless arrays but must be **consistent** (e.g., LJ with σ=1).

## Periodic boundaries

- Use `wrap_positions!` after moving particles.
- The minimum-image convention is implemented by `minimum_image!`.

## Neighbor list rebuild policy

`maybe_rebuild!` triggers a rebuild when the **maximum displacement** since the
last build exceeds `skin/2`. With the **cell-linked builder**, rebuilds are O(N).

## Box size and geometry

Ensure `L > 2*(cutoff + skin)` to avoid ambiguous minimum-image shells.
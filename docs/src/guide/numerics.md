# Numerics & Pitfalls  {#numerics}

## Timestep stability

- Velocity Verlet is symplectic and stable for sufficiently small `dt`.
- If you see exploding energies or trajectories, reduce `dt`.

## Precision

Use `Float64` across positions, velocities, and masses to limit drift.

## Shapes & broadcasting

- Arrays are `(N×D)` with rows as particles.
- `masses` is length `N`. Broadcasting `forces(r) ./ masses` divides each row by its particle mass.

## Boundary conditions (not included)

No periodic or reflective boundaries are included. Implement them in your force or position update logic if needed.

## Potential energy protocol

- Provide `forces(r; return_potential=true) => (F, U)` to enable `potential_energy`.
- If unknown, prefer returning only `F` and skip potential-energy reporting.
# Spec: Module `Verlet.Electrostatics`

Purpose: Placeholder for long-range electrostatics algorithms (Ewald/PME/Wolf).

## Status

- Module exists and includes `ewald.jl`, which currently serves as a scaffold.
- No public API is exported yet; future work may add parameter types and force kernels integrated with `Neighbors` and `Core.ForceField`.

## Planned APIs (non-binding)

- `EwaldReal`, `EwaldReciprocal` parameter types and a composite `EwaldSum` potential.
- PME FFT-based reciprocal-space evaluator hooks.
- Compatibility with `PairTable`/type-pair parameterization.


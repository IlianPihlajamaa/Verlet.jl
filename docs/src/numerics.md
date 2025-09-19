# Numerical Notes

## Constraint Tolerances and Thermostat Parameters

When using [`LangevinBAOABConstrained`](@ref Verlet.Integrators.LangevinBAOABConstrained), the SHAKE and RATTLE solvers enforce
distance constraints up to the tolerance specified in [`DistanceConstraints`](@ref).

### Interplay with `γ` and `T`

* Larger friction `γ` and higher temperature `T` produce stronger stochastic kicks in the OU step.
* This may increase the number of iterations needed for constraint convergence.
* If the solver fails to converge within `maxiter`, consider relaxing `tol` slightly.

### Mass Variance

Constraints involving very light particles (e.g. hydrogens) can converge more slowly
if highly coupled with heavier atoms. This is expected and usually acceptable for
simple bond constraints.

### Degrees of Freedom and Temperature

The function [`instantaneous_temperature`](@ref) uses the full degrees of freedom by default.
For constrained systems, users who want diagnostics consistent with statistical mechanics
should instead call:

```julia
degrees_of_freedom(sys; constraints=cons, remove_com=false)
```

and use this value in temperature calculations.

### Reproducibility

The stochastic step uses the RNG stored in the integrator instance.
Set a fixed seed (e.g. `rng = MersenneTwister(1234)`) in tests or benchmarks for reproducibility.

---

See also: [`constraint_residuals`](@ref) for monitoring how well constraints are satisfied.

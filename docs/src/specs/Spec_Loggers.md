# Spec: Module `Verlet.Loggers`

Purpose: Provide reusable logging infrastructure for molecular simulations. Loggers collect values from `Observable` objects (defined in `Verlet.Core`) into preallocated buffers and manage time and ensemble statistics such as static and time correlations.

## Observables

Observables live in `Verlet.Core`:

- `abstract type Observable end`
- `observe!(obs::Observable, system::System, step::Integer)`
  - Implementations compute the value at `step` and return it (or a mutable buffer reused across calls).
  - Observables never own the long-term storage for their data; loggers are responsible for persisting samples.
- `out_type(obs)::DataType` returns the type of an observation and is used to type-stabilise logger buffers.
- Optional helper: `observed_quantity(obs)::Symbol` → canonical name used by `Loggers` for indexing.

## Logger hierarchy

- `abstract type AbstractLogger end`
- Common API:
  - `step!(logger, system, step)` records data for `step`, usually by delegating to each observable if `step in when_to_log(logger)`.
  - `data(logger)` returns the saved data. Result may depend on the type of logger.
  - `save(file, logger)` saves the data to text file.
  
### `ObservableLogger`

- Stores a single observable's raw time series in logger-owned storage.
- Configuration: `ObservableLogger(observable::Observable; when_to_log = Int[], stride = 1, max_samples = nothing)`.
  - `max_samples` limits the retained history to the most recent `max_samples` entries (useful for rolling buffers).
- Workflow:
  1. `step!` checks whether the current step matches the stride or any explicit `when_to_log` entry.
  2. When sampling, `observe!` is invoked and the returned value is appended to the logger's internal buffer.
  3. Step indices are tracked by the logger, deriving timing from the stride and `when_to_log` schedule.
- Example usage (pseudocode, not implemented):

```julia
obs = KineticEnergy()
logger = ObservableLogger(obs; stride = 1000)
for step in 1:10_000
    integrate_step!(system)
    step!(logger, system, step)
end
```

- Targets equal-time correlation matrices between multiple observables.
- Maintains running totals for `ΣA_i` and `Σ(A_i†A_j)` so that the means and correlations can be materialised on demand without keeping every snapshot.
- Configuration: `StaticCorrelationLogger((obsA, obsB, ...); stride=1, when_to_log=[], product=default)`.
  - `product` is the bilinear map used for the correlation (defaults to multiplication for scalars or `dot` for vector-like data). Correlation matrices retain the element type implied by this map (e.g. complex/dual/unitful values).
- API additions:
  - `means(logger)::NTuple` provides the means `⟨A_1⟩`, `⟨A_2⟩`, etc. (elements are `nothing` if no samples were recorded).
  - `covars(logger)::Matrix` provides the covariances with the promoted element type coming from `product`.
- Designed for large ensembles where storing full traces is too costly.
- Example (pseudocode):

If the observables return vectors (e.g. a value for every particle), an inner product is implied.

```julia
obsA = ParticleVelocity()
obsB = ParticleForce()
logger = StaticCorrelationLogger((obsA, obsB)) # specify in tuple for type stability
logger = StaticCorrelationLogger(obsA) # autocorrelation function

for step in 1:5_000
    step!(integrator, system)
    step!(logger, system, step)
end
```

### `TimeCorrelationLogger`

- Implements the multiple-windows/order-`O(N)` algorithm for auto- or cross-time correlations.
- Configuration: `TimeCorrelationLogger(obsA, obsB; block_length=32, nblocks=6, min_stride=1, timestep=Δt, product=default)`.
  - `obsB` is optional; omitting it computes autocorrelations.
  - `block_length` sets the number of samples stored per window, `nblocks` controls how many stride decades are tracked, and `min_stride` defines the smallest sampling interval.
- `product` is the bilinear map used for the correlation. By default it multiplies scalars and calls `dot` for vector-like data; pass a custom function when working with tensors, units, or AD dual numbers to preserve types.
- Each block keeps ring buffers for the last `block_length` samples at stride `min_stride * block_length^(b-1)` and accumulates sums/counts for every admissible lag so that the full correlation curve is recovered without quadratic work.
- `data(logger)` returns `(; times, values, counts, block_strides, block_length, sample_stride)` where
  - `times` are in physical units (`lag * stride * timestep`),
  - `values` keep the correlation element type inferred from the observables and `product`, and
  - `counts` record how many samples contributed to each lag (useful for error estimates).
- Example (pseudocode):

```julia
vel = VelocityObservable()
logger = TimeCorrelationLogger(vel; block_length=16, nblocks=5, min_stride=10, timestep=Δt)
for step in 0:nsteps
    integrate_step!(system)
    step!(logger, system, step)
end
acf = data(logger)
plot(acf.times, acf.values)
```

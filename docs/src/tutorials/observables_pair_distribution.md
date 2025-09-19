# Tutorial · Observables & Loggers

This tutorial demonstrates the new logging infrastructure. We will implement a simple radial distribution function (RDF) observable and record it with `ObservableLogger` while running a short trajectory.

The design splits responsibilities:

- An `Observable` computes a value for the current system state and returns it. It does **not** manage persistent storage.
- A logger (here `ObservableLogger`) owns the sample buffers, schedules when measurements happen, and stores step indices and values.

## Implementing an RDF observable

An RDF measures pair separations in a histogram. The observable keeps the histogram edges and a reusable workspace that is reset every time it is evaluated.

```julia
using Verlet, StaticArrays

struct RDFObservable <: Observable
    edges::Vector{Float64}
    counts::Vector{Float64}
end

function RDFObservable(rmax::Float64, nbins::Int)
    edges = collect(range(0, rmax; length = nbins + 1))
    return RDFObservable(edges, zeros(Float64, nbins))
end

function observe!(obs::RDFObservable, system::System, step::Integer)
    fill!(obs.counts, 0)
    positions = system.positions
    box = system.box
    n = natoms(system)

    for i in 1:n-1
        ri = positions[i]
        for j in i+1:n
            rj = positions[j]
            Δ = minimum_image(ri - rj, box)
            dist = sqrt(sum(abs2, Δ))
            bin = searchsortedlast(obs.edges, dist)
            if 1 <= bin <= length(obs.counts)
                obs.counts[bin] += 2 # symmetry factor
            end
        end
    end

    return obs.counts
end

out_type(::RDFObservable) = Vector{Float64}
observed_quantity(::RDFObservable) = :g_r
```

`observe!` fills the counts array and returns it. The logger copies the vector before storing it, so the observable stays reusable.

## Running a simulation with logging

```julia
box = CubicBox(5.0)
pos = [SVector{3,Float64}(randn(3)) for _ in 1:32]
vel = [SVector{3,Float64}(0.0, 0.0, 0.0) for _ in 1:32]
forces = [SVector{3,Float64}(0.0, 0.0, 0.0) for _ in 1:32]
masses = ones(Float64, 32)
types = ones(Int, 32)
type_names = Dict(1 => :A)
ff = ForceField(())

system = System(pos, vel, forces, masses, box, types, type_names; forcefield = ff)

integrator = VelocityVerlet(1e-3)
rdf = RDFObservable(2.5, 50)
logger = ObservableLogger(rdf; stride = 10)

for step in 1:1_000
    step!(integrator, system)
    step!(logger, system, step)
end

samples = data(logger)
@show length(samples.steps)        # number of stored frames
@show samples.steps[1:3]
@show samples.values[1]
```

The resulting histogram vector can be normalised after the run. The logger reports which steps were recorded; you can use that metadata to match other time series.

## Combining with other observables

You can simultaneously record built-in observables such as `DensityObservable()` or `TemperatureObservable()` in additional loggers. For composite statistics, use `StaticCorrelationLogger((KineticEnergyObservable(), DensityObservable()))` to keep running means and cross-correlations without storing entire trajectories.

To analyse dynamical quantities, switch to `TimeCorrelationLogger`. Example: the velocity autocorrelation function (VACF)

```julia
vacf = TimeCorrelationLogger(VelocityObservable(); block_length = 32, nblocks = 6,
                              min_stride = 5, timestep = Δt)
for step in 0:nsteps
    step!(integrator, system)
    step!(vacf, system, step)
end

vacf_data = data(vacf)
plot(vacf_data.times, vacf_data.values)
```

The logger returns the non-uniform time grid, averaged correlation values (keeping the same units/AD types as the inputs), and the number of samples per lag so you can weight error bars or discard under-sampled tails. Supply a custom `product = (a, b) -> ...` when you need tensor contractions beyond the default scalar/vector handling.

## Saving results

Both logger types implement `save(file, logger)`, which writes a tab-separated text file. For the RDF logger above:

```julia
save("rdf.txt", logger)
```

Each line contains the step index followed by the histogram entries.

## Next steps

- Add more sophisticated buffering (e.g. block averaging) by implementing custom logger types.
- Use different schedules with `when_to_log=[50, 100, 150]` to capture infrequent events without tightening the stride.

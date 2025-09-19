"""
    abstract type AbstractIntegrator end

Marker supertype for all integration / minimisation schemes.
"""
abstract type AbstractIntegrator end


"""
    integrate!(integrator::AbstractIntegrator, system::System, nsteps::Integer, args...;
               callback=nothing, neighbor_kwargs...) -> System

Run `integrator` against `system` for `nsteps` iterations. Concrete integrators
must implement this method; the default definition throws. When provided, the
`callback` is invoked after each `step!` with `(system, step, integrator)` and
may return `false` to request early termination.
"""
function integrate!(integrator::AbstractIntegrator, system::System, nsteps::Integer, args...;
                    callback=nothing, neighbor_kwargs...)
    nsteps < 0 && throw(ArgumentError("nsteps must be non-negative"))
    nsteps == 0 && return system

    for step in 1:nsteps
        step!(integrator, system)

        if !isnothing(callback)
            cb_result = callback(system, step, integrator)
            cb_result === false && break
        end
        stop_requested(integrator) && break
    end

    return system
end

function step!(integrator::AbstractIntegrator, system::System)
    throw(MethodError(step!, (integrator, system)))
end


function rebuild_neighbors! end

stop_requested(::AbstractIntegrator) = false

function maybe_rebuild(system::System, args...; kwargs...)
    rebuild_neighbors!(system, args...; kwargs...)
    return args[end]
end

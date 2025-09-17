"""
    abstract type AbstractIntegrator end

Marker supertype for all integration / minimisation schemes.
"""
abstract type AbstractIntegrator end

"""
    integrate!(integrator::AbstractIntegrator, system::System, nsteps::Integer, args...; kwargs...) -> System

Run `integrator` against `system` for `nsteps` iterations. Concrete integrators
must implement this method; the default definition throws.
"""
function integrate!(integrator::AbstractIntegrator, system::System, nsteps::Integer, args...; kwargs...)
    throw(MethodError(integrate!, (integrator, system, nsteps, args...); kwargs...))
end

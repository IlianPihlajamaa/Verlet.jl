"""
    VelocityVerlet(dt; wrap=false)

Velocity-Verlet integrator that advances a `System` in time using the
`ForceField` attached to the system. Set `wrap=true` to apply periodic wrapping
after each position update. The optional `callback(system, step)` may return
`false` to terminate the integration early.
"""
struct VelocityVerlet{T} <: AbstractIntegrator
    dt::T
    wrap::Bool
end

VelocityVerlet(dt::Real; wrap::Bool=false) = VelocityVerlet{typeof(dt)}(dt, wrap)

function integrate!(integrator::VelocityVerlet, system::System, nsteps::Integer, args...; callback=nothing)
    isempty(args) || throw(ArgumentError("VelocityVerlet does not accept extra positional arguments"))
    nsteps < 0 && throw(ArgumentError("nsteps must be non-negative"))
    nsteps == 0 && return system

    dt = integrator.dt
    wrap = integrator.wrap

    positions = system.positions
    velocities = system.velocities
    masses = system.masses

    _update_forces!(system)

    for step in 1:nsteps
        current_forces = system.forces
        @inbounds for i in eachindex(positions)
            accel = current_forces[i] / masses[i]
            velocities[i] += 0.5 * dt * accel
            positions[i] += dt * velocities[i]
        end
        wrap && wrap_positions!(positions, system.box)

        _update_forces!(system)
        new_forces = system.forces
        @inbounds for i in eachindex(positions)
            accel = new_forces[i] / masses[i]
            velocities[i] += 0.5 * dt * accel
        end

        if callback !== nothing && callback(system, step) === false
            break
        end
    end

    return system
end

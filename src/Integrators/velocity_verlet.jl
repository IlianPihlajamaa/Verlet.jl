"""
    VelocityVerlet(dt; wrap=false)

Velocity-Verlet integrator that advances a `System` in time using the
`ForceField` attached to the system. Set `wrap=true` to apply periodic wrapping
after each position update. The optional callback supplied to `integrate!`
receives `(system, step, integrator)` and may return `false` to terminate early.
"""
struct VelocityVerlet{T} <: AbstractIntegrator
    dt::T
    wrap::Bool
end

VelocityVerlet(dt::Real; wrap::Bool=false) = VelocityVerlet{typeof(dt)}(dt, wrap)

function step!(integrator::VelocityVerlet, system::System)
    dt = integrator.dt
    wrap = integrator.wrap

    positions = system.positions
    velocities = system.velocities
    masses = system.masses

    forces = _update_forces!(system)

    @inbounds for i in eachindex(positions)
        accel = forces[i] / masses[i]
        velocities[i] += 0.5 * dt * accel
        positions[i] += dt * velocities[i]
    end
    wrap && wrap_positions!(positions, system.box)

    new_forces = _update_forces!(system)
    @inbounds for i in eachindex(positions)
        accel = new_forces[i] / masses[i]
        velocities[i] += 0.5 * dt * accel
    end

    return true
end

# Implementation for integrators, e.g. velocity_verlet!

"""
    velocity_verlet!(system::System{T,IT,Dims}, forces::Function, dt::T) where {T,IT,Dims}

Advance `system` by one timestep `dt` using the velocity Verlet integrator.
"""
function velocity_verlet!(system::System{T,IT,Dims}, forces::Function, dt::T) where {T,IT,Dims}
    dt2 = dt * dt
    current_forces = forces(system.positions)
    accel = current_forces ./ system.masses
    system.positions .+= system.velocities .* dt .+ 0.5 .* accel .* dt2
    new_forces = forces(system.positions)
    new_accel = new_forces ./ system.masses
    system.velocities .+= 0.5 .* (accel .+ new_accel) .* dt
    return system
end

"""
    integrate!(integrator!, system::System, forces, dt, nsteps[, args...]; callback=nothing, kwargs...)

Repeatedly apply `integrator!` to `system` for `nsteps` steps.

- `integrator!` is invoked as `integrator!(system, forces, dt, args...; kwargs...)`.
- `callback`, if provided, is called after each step with `(system, step)`; returning
  `false` stops the loop early.
- `nsteps` must be a non-negative integer. When `nsteps == 0`, the system is returned
  unchanged and the callback is not invoked.
"""
function integrate!(integrator!, system::System, forces, dt, nsteps::Integer, args...; callback=nothing, kwargs...)
    nsteps < 0 && throw(ArgumentError("nsteps must be non-negative, got $nsteps"))
    nsteps == 0 && return system

    if callback === nothing
        for _ in 1:nsteps
            integrator!(system, forces, dt, args...; kwargs...)
        end
    else
        for step in 1:nsteps
            integrator!(system, forces, dt, args...; kwargs...)
            callback(system, step) === false && break
        end
    end

    return system
end

"""
    potential_energy(system::System{T,IT,Dims}, forces::Function) -> T

Try to obtain total potential energy using `forces(...; return_potential=true)` convention.
"""
function potential_energy(system::System{T,IT,Dims}, forces::Function) where {T,IT,Dims}
    F_U_kw = try
        forces(system.positions; return_potential=true)
    catch e
        nothing
    end
    if F_U_kw !== nothing && F_U_kw isa Tuple && length(F_U_kw) == 2
        _, U = F_U_kw
        return T(U)
    end
    error("Force function does not support `return_potential=true`; cannot compute potential_energy.")
end

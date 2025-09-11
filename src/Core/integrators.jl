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
    potential_energy(system::System{T,IT,Dims}, forces::Function) -> T

Try to obtain total potential energy using `forces(...; return_potential=true)` convention.
"""
function potential_energy(system::System{T,IT,Dims}, forces::Function) where {T,IT,Dims}
    F_U_kw = try
        forces(system.positions; return_potential=true)
    catch e
        println("Caught error in potential_energy: $e")
        nothing
    end
    if F_U_kw !== nothing && F_U_kw isa Tuple && length(F_U_kw) == 2
        _, U = F_U_kw
        return T(U)
    end
    error("Force function does not support `return_potential=true`; cannot compute potential_energy.")
end

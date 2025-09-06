# Implementation for integrators, e.g. velocity_verlet!

"""
    velocity_verlet!(system::ParticleSystem, forces::Function, dt::Float64)

Advance `system` by one timestep `dt` using the velocity Verlet integrator.
"""
function velocity_verlet!(system::ParticleSystem, forces::Function, dt::Float64)
    @assert size(system.positions) == size(system.velocities) "positions/velocities must be same size"
    @assert length(system.masses) == size(system.positions, 1) "length(masses) must equal number of particles"
    dt2 = dt * dt
    a = forces(system.positions) ./ system.masses
    system.positions .+= system.velocities .* dt .+ 0.5 .* a .* dt2
    a_new = forces(system.positions) ./ system.masses
    system.velocities .+= 0.5 .* (a .+ a_new) .* dt
    return system
end

"""
    potential_energy(system::ParticleSystem, forces::Function) -> Float64

Try to obtain total potential energy using `forces(...; return_potential=true)` convention.
"""
function potential_energy(system::ParticleSystem, forces::Function)::Float64
    F_U_kw = try
        forces(system.positions; return_potential=true)
    catch
        nothing
    end
    if F_U_kw !== nothing && F_U_kw isa Tuple && length(F_U_kw) == 2
        _, U = F_U_kw
        return float(U)
    end
    error("Force function does not support `return_potential=true`; cannot compute potential_energy.")
end

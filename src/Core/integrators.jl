# Implementation for integrators, e.g. velocity_verlet!

"""
    velocity_verlet!(system::ParticleSystem{Dims,T_float}, forces::Function, dt::T_float)

Advance `system` by one timestep `dt` using the velocity Verlet integrator.
"""
function velocity_verlet!(system::ParticleSystem{Dims,T_float}, forces::Function, dt::T_float) where {Dims,T_float}
    @assert length(system.positions) == length(system.velocities) "positions/velocities must be same size"
    @assert length(system.masses) == length(system.positions) "length(masses) must equal number of particles"
    dt2 = dt * dt
    a = forces(system.positions) 
    a ./= system.masses
    system.positions .= system.positions .+ system.velocities .* dt .+ a .* (0.5 * dt2)
    a_new = forces(system.positions)
    a./= system.masses
    system.velocities .= system.velocities .+ 0.5 .* (a .+ a_new) .* dt
    return system
end

"""
    potential_energy(system::ParticleSystem{Dims,T_float}, forces::Function) -> T_float

Try to obtain total potential energy using `forces(...; return_potential=true)` convention.
"""
function potential_energy(system::ParticleSystem{Dims,T_float}, forces::Function) where {Dims,T_float}
    F_U_kw = try
        forces(system.positions; return_potential=true)
    catch e 
        println("Caught error in potential_energy: $e")
        nothing
    end
    if F_U_kw !== nothing && F_U_kw isa Tuple && length(F_U_kw) == 2
        _, U = F_U_kw
        return T_float(U)
    end
    error("Force function does not support `return_potential=true`; cannot compute potential_energy.")
end



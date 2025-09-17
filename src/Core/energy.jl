"""
    potential_energy(system::System{T}, forces::Function) -> T

Obtain the potential energy by calling `forces(system.positions; return_potential=true)`.
If the force function does not support this protocol, the helper will raise an
error.
"""
function potential_energy(system::System{T,IT,D}, forces::Function) where {T,IT,D}
    F_U_kw = try
        forces(system.positions; return_potential=true)
    catch
        nothing
    end
    if F_U_kw !== nothing && F_U_kw isa Tuple && length(F_U_kw) == 2
        _, U = F_U_kw
        return T(U)
    end
    error("Force function does not support `return_potential=true`; cannot compute potential_energy.")
end

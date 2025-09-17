function _ensure_forcefield(system::System)
    ff = system.forcefield
    ff === nothing && throw(ArgumentError("system has no forcefield assigned"))
    return ff
end

function _update_forces!(system::System)
    ff = _ensure_forcefield(system)
    compute_all_forces!(system, ff)
    return system.forces
end

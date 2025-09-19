function compute_forces!(pot::AbstractPotential, sys::System)
    throw(MethodError(compute_forces!, (pot, sys)))
end

function prepare_neighbors!(ff, sys::System)
    return ff
end

function compute_all_forces!(sys::System, ff)
    fill!(sys.forces, zero(eltype(sys.forces)))
    prepare_neighbors!(ff, sys)

    for pot in ff.layers
        compute_forces!(pot, sys)
    end

    for interaction in sys.specific_potentials
        compute_forces!(interaction, sys)
    end
    return sys
end

function compute_all_forces!(sys::System)
    ff = sys.forcefield
    ff === nothing && throw(ArgumentError("System has no forcefield assigned"))
    compute_all_forces!(sys, ff)
end

function compute_potential_energy(sys::System, ff)
    E = zero(eltype(sys.positions[1]))
    prepare_neighbors!(ff, sys)

    for pot in ff.layers
        E += compute_potential_energy(pot, sys)
    end

    for interaction in sys.specific_potentials
        E += compute_potential_energy(interaction, sys)
    end

    return E
end

function compute_potential_energy(sys::System)
    ff = sys.forcefield
    ff === nothing && throw(ArgumentError("System has no forcefield assigned"))
    compute_potential_energy(sys, ff)
end

function compute_potential_energy(pot::AbstractPotential, sys::System)
    return zero(eltype(sys.positions[1]))
end

function compute_forces! end
function compute_all_forces! end

struct ForceField{ForcesTuple}
    layers::ForcesTuple
end

function compute_all_forces!(sys::System, ff::ForceField)
    # Reset forces before calculation
    fill!(sys.forces, zero(eltype(sys.forces)))

    # Pairwise forces from the force field layers
    for pot in ff.layers
        compute_forces!(pot, sys)
    end

    # Specific bonded forces
    for interaction in sys.specific_potentials
        compute_forces!(interaction, sys)
    end
end

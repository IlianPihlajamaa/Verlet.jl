using StructArrays

struct ForceField{ForcesTuple}
    layers::ForcesTuple
end

function is_excluded(pot::AbstractPairPotential, i::Integer, j::Integer)
    # This is a placeholder. A real implementation would be more efficient.
    return (i, j) in pot.exclusions || (j, i) in pot.exclusions
end

function build_neighbors_from_master!(pot::AbstractPairPotential, sys::System, masterNL::MasterNeighborList)
    PairType = eltype(pot.params.table)

    neighbors = pot.neighbors
    # empty! on a StructArray clears the underlying vectors, but preserves their capacity.
    # This allows for non-allocating rebuilds if the capacity is sufficient.
    empty!(neighbors)

    sizehint!(neighbors, length(masterNL.entries))

    for entry in masterNL.entries
        p = pot.params.table[sys.types[entry.i], sys.types[entry.j]]
        if !is_excluded(pot, entry.i, entry.j) && entry.r2 < (p.rc + pot.skin)^2
            push!(neighbors, (i=entry.i, j=entry.j, pair=p))
        end
    end
end


function build_all_neighbors!(master_nl::MasterNeighborList, ff::ForceField, sys::System; method::Symbol=:cells)
    # 1. Compute master NL radius
    rc_max = maximum(
        maximum([p.rc for p in pot.params.table]) + pot.skin
        for pot in ff.layers
    )

    build_master_neighborlist!(master_nl, sys.positions, sys.box; r_verlet=rc_max, method=method)

    # 2. Build per-potential NLs
    for pot in ff.layers
        build_neighbors_from_master!(pot, sys, master_nl)
    end
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

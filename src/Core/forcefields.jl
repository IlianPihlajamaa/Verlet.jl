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

    ivec = T_int[]
    jvec = T_int[]
    pairvec = PairType[]
    sizehint!(ivec, length(masterNL.entries))
    sizehint!(jvec, length(masterNL.entries))
    sizehint!(pairvec, length(masterNL.entries))

    for entry in masterNL.entries
        p = pot.params.table[sys.types[entry.i], sys.types[entry.j]]
        if !is_excluded(pot, entry.i, entry.j) && entry.r2 < (p.rc + pot.skin)^2
            push!(ivec, entry.i)
            push!(jvec, entry.j)
            push!(pairvec, p)
        end
    end
    pot.neighbors = StructArray{NeighborPair{PairType, T_int}}((i=ivec, j=jvec, pair=pairvec))
end


function build_all_neighbors!(ff::ForceField, sys::System, master_skin::T_Float; method::Symbol=:cells)
    # 1. Compute master NL radius
    rc_max = maximum(
        maximum([p.rc for p in pot.params.table]) + pot.skin
        for pot in ff.layers
    )

    masterNL = build_master_neighborlist(sys.positions, sys.box; r_verlet=rc_max, skin=master_skin, method=method)

    # 2. Build per-potential NLs
    for pot in ff.layers
        build_neighbors_from_master!(pot, sys, masterNL)
    end

    return masterNL
end

function compute_all_forces!(sys::System, ff::ForceField)
    # Reset forces before calculation
    fill!(sys.forces, zero(eltype(sys.forces)))

    for pot in ff.layers
        compute_forces!(pot, sys)
    end
end

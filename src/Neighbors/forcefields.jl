
function is_excluded(pot::AbstractPairPotential, i::Integer, j::Integer)
    # This is a placeholder. A real implementation would be more efficient.
    return (i, j) in pot.exclusions
end

function build_neighbors_from_master!(pot::AbstractPairPotential, sys::System, masterNL::MasterNeighborList)
    PairType = eltype(pot.params.table)

    neighbors = pot.neighborlist.neighbors
    empty!(neighbors)

    types = sys.types
    skin = pot.skin
    positions = sys.positions
    box = sys.box
    for entry in masterNL.pairs
        i, j = entry[1], entry[2]
        ri = positions[i]
        rj = positions[j]
        Δ = ri - rj
        Δ = minimum_image(Δ, box)
        r2 = dot(Δ, Δ)
        p = pot.params.table[types[i], types[j]]
        rc = p.rc
        if !is_excluded(pot, i, j) && r2 < (rc + skin)^2
            push!(neighbors, NeighborPair{PairType, Int}(i, j, p))
        end
    end
end

function build_all_neighbors!(master_nl::MasterNeighborList, ff::ForceField, sys::System; method::Symbol=:cells)
    layers = ff.layers
    rc_max = maximum(
        maximum(p.rc for p in pot.params.table) + pot.skin
        for pot in layers
    )

    build_master_neighborlist!(master_nl, sys; r_verlet=rc_max, method=method)

    map(pot -> build_neighbors_from_master!(pot, sys, master_nl), layers)
end

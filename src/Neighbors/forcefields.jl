
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
    for entry in masterNL.entries
        i, j, r2 = entry.i, entry.j, entry.r2
        p = pot.params.table[types[i], types[j]]
        rc = p.rc
        if !is_excluded(pot, i, j) && r2 < (rc + skin)^2
            push!(neighbors, NeighborPair{PairType, Int}(i, j, p))
        end
    end
end



function build_neighbors_from_master!(pot::AbstractPairPotential, sys::System, masterNL::MasterNeighborCSRList)
    PairType = eltype(pot.params.table)

    neighbors = pot.neighborlist.neighbors
    empty!(neighbors)

    types = sys.types
    skin = pot.skin
    R = sys.positions
    box = sys.box

    rowptr = masterNL.rowptr
    colind = masterNL.colind
    N = length(R)

    @inbounds for i in 1:N
        starti = rowptr[i]
        endi = rowptr[i+1] - 1
        for k in starti:endi
            j = colind[k]
            p = pot.params.table[types[i], types[j]]
            rc = p.rc
            Δ = R[i] - R[j]
            Δ = minimum_image(Δ, box)
            r2 = dot(Δ, Δ)
            if !is_excluded(pot, i, j) && r2 < (rc + skin)^2
                push!(neighbors, NeighborPair{PairType, Int}(i, j, p))
            end
        end
    end
end

function build_all_neighbors!(master_nl::MasterNeighborList, ff::ForceField, sys::System; method::Symbol=:cells)
    # 1. Compute master NL radius
    layers = ff.layers

    rc_max = maximum(
        maximum([p.rc for p in pot.params.table]) + pot.skin
        for pot in layers
    )

    build_master_neighborlist!(master_nl, sys.positions, sys.box; r_verlet=rc_max, method=method)

    # 2. Build per-potential NLs
    map(pot -> build_neighbors_from_master!(pot, sys, master_nl), layers)
end


function build_all_neighbors!(master_nl::MasterNeighborCSRList, ff::ForceField, sys::System; method::Symbol=:cells)
    # 1. Compute master NL radius
    layers = ff.layers
    rc_max = maximum(
        maximum([p.rc for p in pot.params.table]) + pot.skin
        for pot in layers
    )

    build_master_neighborlist!(master_nl, sys.positions, sys.box; r_verlet=rc_max, method=method)

    # 2. Build per-potential NLs
    map(pot -> build_neighbors_from_master!(pot, sys, master_nl), layers)
end


function build_all_neighbors!(master_nl::MyNeighborList, ff::ForceField, sys; method::Symbol=:cells)
    # 1. Compute master NL radius
    layers = ff.layers

    rc_max = maximum(
        maximum([p.rc for p in pot.params.table]) + pot.skin
        for pot in layers
    )

    rebuild!(master_nl, sys, sys.box)


    # 2. Build per-potential NLs
    map(pot -> build_neighbors_from_master!(pot, sys, master_nl), layers)
end

function build_neighbors_from_master!(pot::AbstractPairPotential, sys::System, masterNL::MyNeighborList)
    PairType = eltype(pot.params.table)

    neighbors = pot.neighborlist.neighbors
    empty!(neighbors)

    types = sys.types
    skin = pot.skin
    for entry in masterNL.pairs
        i, j = entry[1], entry[2]
        ri = sys.positions[i]
        rj = sys.positions[j]
        Δ = ri - rj
        Δ = minimum_image(Δ, sys.box)
        r2 = dot(Δ, Δ)
        p = pot.params.table[types[i], types[j]]
        rc = p.rc
        if !is_excluded(pot, i, j) && r2 < (rc + skin)^2
            push!(neighbors, NeighborPair{PairType, Int}(i, j, p))
        end
    end
end
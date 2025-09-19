
struct EmptyNeighborList end

mutable struct ForceField{Layers}
    layers::Layers
    master::Any
    master_method::Symbol
end

function ForceField(layers; method::Symbol=:cells)
    master = _initial_master(layers)
    ForceField{typeof(layers)}(layers, master, method)
end

function _initial_master(layers)
    _has_pair_potentials(layers) ? nothing : EmptyNeighborList()
end

function prepare_neighbors!(ff::ForceField, sys::System)
    pair_layers = _collect_pair_potentials(ff.layers)
    isempty(pair_layers) && return ff

    rc_max, skin_max = _pair_extents(pair_layers)
    skin_max = max(skin_max, zero(rc_max))

    master = ff.master
    if master === nothing || master isa EmptyNeighborList
        master = MasterNeighborList(sys; cutoff=rc_max, skin=skin_max)
        ff.master = master
    elseif !(master isa MasterNeighborList)
        master = MasterNeighborList(sys; cutoff=rc_max, skin=skin_max)
        ff.master = master
    end

    force_rebuild = false
    if master isa MasterNeighborList
        if master.skin != skin_max
            master.skin = skin_max
            force_rebuild = true
        end
        if master.cutoff != rc_max
            master.cutoff = rc_max
            force_rebuild = true
        end
    end

    r_verlet = rc_max + skin_max
    if force_rebuild || _needs_master_rebuild(master, sys, skin_max)
        build_master_neighborlist!(master, sys; r_verlet=r_verlet, method=ff.master_method)
    end

    for pot in pair_layers
        build_neighbors_from_master!(pot, sys, master)
    end

    return ff
end

function is_excluded(pot::AbstractPairPotential, i::Integer, j::Integer)
    # This is a placeholder. A real implementation would be more efficient.
    return (i, j) in pot.exclusions
end

function build_neighbors_from_master!(pot::AbstractPairPotential, sys::System, masterNL::MasterNeighborList)
    PairType = eltype(pot.params.table)

    neighbors = pot.neighborlist.neighbors
    empty!(neighbors)
    NeighborType = eltype(neighbors)
    IndexType = fieldtype(NeighborType, :i)

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
            push!(neighbors, NeighborPair{PairType, IndexType}(IndexType(i), IndexType(j), p))
        end
    end
end

function build_all_neighbors!(master_nl::MasterNeighborList, ff::ForceField, sys::System; method::Symbol=:cells)
    pair_layers = [pot for pot in ff.layers if pot isa AbstractPairPotential]
    isempty(pair_layers) && return master_nl

    rc_max = maximum(
        maximum(p.rc for p in pot.params.table) + pot.skin
        for pot in pair_layers
    )

    build_master_neighborlist!(master_nl, sys; r_verlet=rc_max, method=method)

    for pot in pair_layers
        build_neighbors_from_master!(pot, sys, master_nl)
    end

    return master_nl
end

function _has_pair_potentials(layers)
    for pot in layers
        pot isa AbstractPairPotential && return true
    end
    return false
end

function _collect_pair_potentials(layers)
    result = Vector{Any}()
    for pot in layers
        pot isa AbstractPairPotential && push!(result, pot)
    end
    return result
end

function _pair_extents(pots)
    first_pot = first(pots)
    rc_max = zero(first_pot.params.table[1].rc)
    skin_max = zero(first_pot.skin)
    for pot in pots
        rc_max = max(rc_max, _max_cutoff(pot))
        skin_max = max(skin_max, pot.skin)
    end
    return rc_max, skin_max
end

function _max_cutoff(pot)
    max_rc = zero(pot.params.table[1].rc)
    for val in pot.params.table
        max_rc = max(max_rc, val.rc)
    end
    return max_rc
end

function _needs_master_rebuild(master, sys::System, skin_max)
    master isa EmptyNeighborList && return false
    master === nothing && return false

    positions = sys.positions
    if length(master.r0) != length(positions)
        return true
    end

    threshold = skin_max / 2
    threshold2 = threshold > zero(threshold) ? threshold^2 : zero(threshold)

    max_disp2 = zero(eltype(master.r0[1][1]))
    @inbounds for i in eachindex(positions)
        disp = positions[i] - master.r0[i]
        disp = minimum_image(disp, sys.box)
        max_disp2 = max(max_disp2, dot(disp, disp))
        if max_disp2 > threshold2
            master.max_disp2 = max(master.max_disp2, max_disp2)
            return true
        end
    end
    master.max_disp2 = max(master.max_disp2, max_disp2)
    return false
end

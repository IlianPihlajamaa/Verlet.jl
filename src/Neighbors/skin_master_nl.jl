
using LinearAlgebra
using StaticArrays
import ..Core
import ..Core: System, rebuild_neighbors!, maybe_rebuild


@inline function _wrap_frac!(sf::MVector{D,T}) where {D,T}
    @inbounds for d in 1:D
        x = sf[d]
        sf[d] = x - floor(x)
    end
    return sf
end

@inline _clamp_cell_index(idx::Int, n::Int) = (idx < 0 ? 0 : (idx >= n ? n - 1 : idx))

@inline function _wrap_cell_index(idx::Int, n::Int)
    return idx >= 0 ? (idx % n) : ((n - (-idx % n)) % n)
end

# 0-based → 1-based linear index
@inline function _linid0(idx::NTuple{D,Int}, dims::NTuple{D,Int}) where {D}
    id = 1
    stride = 1
    @inbounds for d in 1:D
        id += (idx[d]) * stride
        stride *= dims[d]
    end
    return id
end

@generated function _neighbor_offsets(::Val{D}) where {D}
    offsets = NTuple{D,Int}[]
    current = Vector{Int}(undef, D)
    function rec(depth)
        if depth > D
            push!(offsets, tuple(current...))
        else
            for x in (-1, 0, 1)
                current[depth] = x
                rec(depth + 1)
            end
        end
        return nothing
    end
    rec(1)
    exprs = map(offsets) do tup
        Expr(:tuple, tup...)
    end
    return Expr(:tuple, exprs...)
end
"""
    MasterNeighborList{D,T}

Master Verlet neighbor list storing unique `(i, j)` candidates within
`cutoff + skin` for a `D`-dimensional periodic system. The list retains
buffers used by the cell-linked builder so repeated rebuilds avoid
allocations.

Create instances via [`MasterNeighborList(sys; cutoff, skin)`](@ref) or
[`MasterNeighborList(positions, box; cutoff, skin)`](@ref).
"""
mutable struct MasterNeighborList{D,T}
    cutoff::T
    cutoff2::T
    skin::T
    rskin2::T
    pairs::Vector{SVector{2,Int}}   # candidate (i<j) pairs within cutoff+skin
    r0::Vector{SVector{D,T}}        # reference positions at last rebuild
    max_disp2::T
    nbuilds::Int
    ncells::NTuple{D,Int}
    # reusable buffers for cell lists
    head::Vector{Int}
    next::Vector{Int}
end


function Base.show(io::IO, nl:: MasterNeighborList{D,T}) where {D,T}
    print(io,
        " MasterNeighborList{$D,$(T)}(",
        "pairs=", length(nl.pairs),
        ", cutoff=", nl.cutoff,
        ", skin=", nl.skin,
        ", nbuilds=", nl.nbuilds,
        ", cells=", nl.ncells,
        ")")
end

# --- Cell grid choice -------------------------------------------------

@inline function _choose_cell_grid(box, rskin::T, ::Val{3}) where T
    L = box.L
    n = max(1, Int(floor(L / rskin)))
    return (n, n, n)
end

@inline function _choose_cell_grid(box, rskin::T, ::Val{D}) where {T,D}
    L = box.L
    n = max(1, Int(floor(L / rskin)))
    return ntuple(_ -> n, D)
end

function MasterNeighborList(sys::System; cutoff::T, skin::T) where {T<:AbstractFloat}
    _make_neighbor_list(sys.positions, sys.box; cutoff, skin)
end

function MasterNeighborList(positions, box; cutoff::T, skin::T) where {T<:AbstractFloat}
    _make_neighbor_list(positions, box; cutoff, skin)
end

# Main constructor specialized on Box{...,ORTHO}
function _make_neighbor_list(
    positions,
    box;
    cutoff::T, skin::T
) where {T}
    cutoff2 = cutoff * cutoff
    rskin  = cutoff + skin
    rskin2 = rskin * rskin

    N  = length(positions)
    D = length(positions[1])
    r0 = Vector{SVector{D,T}}(undef, N)
    @inbounds for i in 1:N
        r0[i] = positions[i]
    end

    dims = _choose_cell_grid(box, rskin, Val(D))
    ncelltot = prod(dims)

    head = fill(0, max(ncelltot, 1))
    next = Vector{Int}(undef, N)

    pairs = Vector{SVector{2,Int}}()
    sizehint!(pairs, max(16, 100N))

    _build_pairs_cells!(pairs, head, next, positions, box, rskin, dims, Val(D))


    return  MasterNeighborList{D,T}(
        cutoff, cutoff2, skin, rskin2,
        pairs, r0, zero(T), 1, dims, 
        head, next
    )
end


function rebuild!(nl::MasterNeighborList{D,T}, sys::System; method::Symbol=:cells, cutoff=nl.cutoff) where {D,T}
    rebuild!(nl, sys.positions, sys.box; method, cutoff)
end

function rebuild!(nl::MasterNeighborList{D,T}, positions, box; method::Symbol=:cells, cutoff=nl.cutoff) where {D,T}
    nl.cutoff = cutoff
    nl.cutoff2 = cutoff * cutoff
    rskin = nl.cutoff + nl.skin
    nl.rskin2 = rskin * rskin

    N = length(positions)
    if length(nl.r0) != N
        resize!(nl.r0, N)
    end

    if method == :cells
        # D = length(nl.ncells)
        ncells = _choose_cell_grid(box, rskin, Val(D))
        nct = max(prod(ncells), 1)
        if length(nl.head) != nct
            resize!(nl.head, nct)
        end
        fill!(nl.head, 0)
        if length(nl.next) != N
            resize!(nl.next, N)
        end
        empty!(nl.pairs)
        nl.ncells = ncells
        _build_pairs_cells!(nl.pairs, nl.head, nl.next, positions, box, rskin, nl.ncells, Val(D))
    elseif method == :bruteforce
        empty!(nl.pairs)
        _build_pairs_bruteforce!(nl.pairs, positions, box, nl.rskin2)
        nl.ncells = ntuple(_ -> 1, length(nl.ncells))
    elseif method == :all_pairs
        empty!(nl.pairs)
        _build_pairs_allpairs!(nl.pairs, positions)
        nl.ncells = ntuple(_ -> 1, length(nl.ncells))
    else
        error("Unknown neighborlist method: $method")
    end

    @inbounds for i in eachindex(positions)
        nl.r0[i] = positions[i]
    end
    nl.max_disp2 = zero(T)
    nl.nbuilds  += 1
    return nl
end

function Core.maybe_rebuild(system::System, nl::MasterNeighborList; kwargs...)
    positions = system.positions
    if length(nl.r0) != length(positions)
        return rebuild_neighbors!(system, nl; kwargs...)
    end

    skin = nl.skin
    skin_half = skin / 2
    threshold2 = skin_half * skin_half

    max_disp2 = _max_displacement2(positions, nl.r0, system.box)
    nl.max_disp2 = max(nl.max_disp2, max_disp2)

    if max_disp2 > threshold2
        return rebuild_neighbors!(system, nl; kwargs...)
    end

    return nl
end

function _max_displacement2(current, reference, box)
    vec_type = eltype(reference)
    scalar_T = eltype(vec_type)
    isempty(current) && return zero(scalar_T)

    max_disp2 = zero(scalar_T)
    @inbounds for i in eachindex(current)
        disp = current[i] - reference[i]
        disp = minimum_image(disp, box)
        max_disp2 = max(max_disp2, dot(disp, disp))
    end
    return max_disp2
end


fractional(box, ri) = box.L \ ri 

function _build_pairs_cells!(
    pairs::Vector{SVector{2,Int}},
    head::Vector{Int},
    next::Vector{Int},
    positions,
    box,
    rskin::T,
    nc::NTuple{3,Int},
    ::Val{3}
) where {T<:AbstractFloat}
    rskin2 = rskin * rskin
    N = length(positions)
    D = 3
    @assert prod(nc) == length(head) "head buffer size mismatch."

    fill!(head, 0)

    @inbounds for i in 1:N
        ri = positions[i]
        sS = fractional(box, ri)
        s = MVector{D,T}(undef)
        @inbounds for d in 1:D
            s[d] = sS[d]
        end
        _wrap_frac!(s)
        idx = ntuple(d -> begin
            sd = s[d] * nc[d]
            k = Int(floor(sd))
            _wrap_cell_index(k, nc[d])
        end, D)
        lid = _linid0(idx, nc)
        next[i] = head[lid]
        head[lid] = i
    end

    empty!(pairs)

    @inbounds for i3 in 0:(nc[3]-1), i2 in 0:(nc[2]-1), i1 in 0:(nc[1]-1)
        lidA = _linid0((i1, i2, i3), nc)
        uniq = Vector{Int}(undef, 27); nuniq = 0
        for d3 in -1:1, d2 in -1:1, d1 in -1:1
            j1 = _wrap_cell_index(i1 + d1, nc[1])
            j2 = _wrap_cell_index(i2 + d2, nc[2])
            j3 = _wrap_cell_index(i3 + d3, nc[3])
            (0 <= j1 < nc[1] && 0 <= j2 < nc[2] && 0 <= j3 < nc[3]) || continue
            lidB = _linid0((j1, j2, j3), nc)
            (lidB >= lidA) || continue
            seen = false
            @inbounds for t in 1:nuniq
                if uniq[t] == lidB
                    seen = true
                    break
                end
            end
            if !seen
                nuniq += 1
                uniq[nuniq] = lidB
            end
        end
        @inbounds for t in 1:nuniq
            _accumulate_cell_pairs!(pairs, head, next, lidA, uniq[t], positions, box, rskin2)
        end
    end

    return pairs
end

function _build_pairs_cells!(
    pairs::Vector{SVector{2,Int}},
    head::Vector{Int},
    next::Vector{Int},
    positions,
    box,
    rskin::T,
    nc::NTuple{D,Int},
    ::Val{D}
) where {D,T<:AbstractFloat}
    rskin2 = rskin * rskin
    N = length(positions)
    @assert prod(nc) == length(head) "head buffer size mismatch."

    fill!(head, 0)

    @inbounds for i in 1:N
        ri = positions[i]
        sS = fractional(box, ri)
        s = MVector{D,T}(undef)
        @inbounds for d in 1:D
            s[d] = sS[d]
        end
        _wrap_frac!(s)
        idx = ntuple(d -> begin
            sd = s[d] * nc[d]
            k = Int(floor(sd))
            _wrap_cell_index(k, nc[d])
        end, D)
        lid = _linid0(idx, nc)
        next[i] = head[lid]
        head[lid] = i
    end

    empty!(pairs)

    dims = ntuple(d -> 0:(nc[d]-1), D)
    offsets = _neighbor_offsets(Val(D))
    uniq = Vector{Int}(undef, length(offsets))

    @inbounds for idx in CartesianIndices(dims)
        tup = Tuple(idx)
        lidA = _linid0(tup, nc)
        nuniq = 0
        for off in offsets
            idxB = ntuple(d -> _wrap_cell_index(tup[d] + off[d], nc[d]), D)
            lidB = _linid0(idxB, nc)
            lidB < lidA && continue
            seen = false
            @inbounds for t in 1:nuniq
                if uniq[t] == lidB
                    seen = true
                    break
                end
            end
            if !seen
                nuniq += 1
                uniq[nuniq] = lidB
            end
        end
        @inbounds for t in 1:nuniq
            _accumulate_cell_pairs!(pairs, head, next, lidA, uniq[t], positions, box, rskin2)
        end
    end

    return pairs
end

# Accumulate i<j pairs from cell A vs cell B, testing against rskin2
@inline function _accumulate_cell_pairs!(pairs, head, next, lidA::Int, lidB::Int, positions, box, rskin2)
    i = head[lidA]
    if lidA == lidB
        @inbounds while i != 0
            j = next[i]; ri = positions[i]
            while j != 0
                r2 = distance2_minimum_image(ri, positions[j], box)
                if r2 <= rskin2
                    push!(pairs, ifelse(i<j, SVector(i, j), SVector(j,i)))
                end
                j = next[j]
            end
            i = next[i]
        end
    else
        @inbounds while i != 0
            j = head[lidB]; ri = positions[i]
            while j != 0
                r2 = distance2_minimum_image(ri, positions[j], box)
                if r2 <= rskin2
                    if i>j
                        push!(pairs, SVector(j,i))
                    else
                        push!(pairs, SVector(i,j))
                    end
                end
                j = next[j]
            end
            i = next[i]
        end
    end
    return nothing
end

@inline function _build_pairs_bruteforce!(pairs::Vector{SVector{2,Int}}, positions, box, rskin2)
    N = length(positions)
    @inbounds for i in 1:max(N - 1, 0)
        ri = positions[i]
        for j in i+1:N
            r2 = distance2_minimum_image(ri, positions[j], box)
            if r2 <= rskin2
                push!(pairs, SVector(i, j))
            end
        end
    end
    return pairs
end

@inline function _build_pairs_allpairs!(pairs::Vector{SVector{2,Int}}, positions)
    N = length(positions)
    @inbounds for i in 1:max(N - 1, 0)
        for j in i+1:N
            push!(pairs, SVector(i, j))
        end
    end
    return pairs
end

function brute_force_pairs(positions, box, cutoff)
    pairs = Vector{SVector{2,Int}}()
    _build_pairs_bruteforce!(pairs, positions, box, cutoff * cutoff)
    return pairs
end

@inline function brute_force_pairs(sys::System, cutoff)
    brute_force_pairs(sys.positions, sys.box, cutoff)
end

"""
    build_master_neighborlist!(nl::MasterNeighborList, sys::System; r_verlet, method=:cells)
    build_master_neighborlist!(nl::MasterNeighborList, positions, box; r_verlet, method=:cells)

Update a [`MasterNeighborList`](@ref) in-place so that it contains all candidate
`(i, j)` pairs consistent with the requested Verlet radius `r_verlet`.

`method` chooses the builder:

* `:cells` (default) uses the cell-linked list sweep and is `O(N)` for large systems.
* `:bruteforce` walks all `i<j` pairs (`O(N^2)`), useful for debugging and validation.
* `:all_pairs` stores every `(i, j)` combination regardless of the cutoff.

The list's `cutoff` is updated to `max(zero(typeof(nl.cutoff)), r_verlet - nl.skin)`
so that the stored `pairs` are within `cutoff + skin`. Returns `nl`.
"""
function build_master_neighborlist!(nl::MasterNeighborList{D,T}, sys::System; r_verlet::Real, method::Symbol=:cells) where {D,T}
    build_master_neighborlist!(nl, sys.positions, sys.box; r_verlet, method)
end

function build_master_neighborlist!(nl::MasterNeighborList{D,T}, positions, box; r_verlet::Real, method::Symbol=:cells) where {D,T}
    if method === :all_pairs
        rebuild!(nl, positions, box; method=:all_pairs, cutoff=convert(T, r_verlet))
    else
        rskin = convert(T, r_verlet)
        cutoff = max(zero(T), rskin - nl.skin)
        rebuild!(nl, positions, box; method=method, cutoff=cutoff)
    end
end

@inline function displacement(ri::SVector{3,T}, rj::SVector{3,T}, box) where {T}
    dx = rj[1] - ri[1]
    dy = rj[2] - ri[2]
    dz = rj[3] - ri[3]
    L = box.L
    dx_corr = round(dx/L)*L
    dy_corr = round(dy/L)*L
    dz_corr = round(dz/L)*L
    SVector{3,T}(dx - dx_corr, dy - dy_corr, dz - dz_corr)
end



@inline function displacement(ri::SVector{D,T}, rj::SVector{D,T}, box) where {D,T}
    dx = rj - ri
    L = box.L
    dx_corr = @inbounds SVector{D,T}(ntuple(d ->  round(dx[d]/L)*L , D))
    dx - dx_corr
end


@inline function distance2_minimum_image(ri::SVector{D,T}, rj::SVector{D,T}, box) where {D,T}
    dr = displacement(ri, rj, box)
    dot(dr, dr)
end

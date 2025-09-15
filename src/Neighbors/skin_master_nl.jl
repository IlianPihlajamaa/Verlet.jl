# =====================================================================
# File: src/Core/ MyNeighborLists.jl
# =====================================================================
# Kinetikos: CPU neighbor lists (cell list + Verlet) with skin & rebuild
# =====================================================================



using LinearAlgebra
using StaticArrays

# Expected aliases:
# const Vec{D,T} = SVector{D,T}
# const Mat{D,T,D2} = SMatrix{D,D,T,D2}
# struct Box{D,T,D2,ORTHO} end
# System must have a field `box`

# --- Utilities --------------------------------------------------------

# Physical lengths of the (possibly skewed) box edges (|a|, |b|, |c|)
@inline function _edge_lengths(h)
    D = size(h, 1)
    if D == 1
        return (abs(h[1, 1]),)
    elseif D == 2
        a = @SVector [h[1, 1], h[2, 1]]
        b = @SVector [h[1, 2], h[2, 2]]
        return (norm(a), norm(b))
    else
        a = @SVector [h[1, 1], h[2, 1], h[3, 1]]
        b = @SVector [h[1, 2], h[2, 2], h[3, 2]]
        c = @SVector [h[1, 3], h[2, 3], h[3, 3]]
        return (norm(a), norm(b), norm(c))
    end
end

# Wrap fractional coordinate to [0,1) for periodic dims
@inline function _wrap_frac!(sf::MVector{D,T}, periodic::Bool) where {D,T}
    @inbounds for d in 1:D
        if periodic
            x = sf[d]
            sf[d] = x - floor(x)
        end
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

# --- Data structure ---------------------------------------------------

"""
     MyNeighborList{D,T}

CPU neighbor/Verlet list for D-dimensional systems (element type `T`).

Fields (hot-path layout; stable across rebuilds):
- `cutoff::T`      : interaction cutoff (used at iteration time)
- `cutoff2::T`     : `cutoff^2` (cached)
- `skin::T`        : Verlet skin (buffer for rebuild policy)
- `rskin2::T`      : `(cutoff + skin)^2` (build radius squared)
- `pairs::Vector{SVector{2,Int}}`
    Candidate (i<j) pairs within `(cutoff + skin)` at last rebuild.
- `r0::Vector{SVector{D,T}}`
    Reference positions at last rebuild (for displacement tracking).
- `max_disp2::T`   : max squared displacement since last rebuild
- `nbuilds::Int`   : number of rebuilds performed (including initial)
- `ncells::NTuple{D,Int}`
    Cell grid shape used at last build (e.g. `(nx,ny,nz)`); `(1,...)` if none.
- `use_cells::Bool`
    Whether the cell builder is desired when the box is orthorhombic.
- `head::Vector{Int}` / `next::Vector{Int}`
    Reusable buffers for the linked-cell list; zero-filled between rebuilds.

Design notes
------------
- For orthorhombic boxes (`ORTHO=true`) and `use_cells=true`, building is
  expected O(N). For triclinic (`ORTHO=false)`, an O(N²) fallback is used.
- Rebuild policy is controlled by `skin`: rebuild iff `2*max_disp > skin`.
"""

mutable struct MyNeighborList{D,T}
    cutoff::T
    cutoff2::T
    skin::T
    rskin2::T
    pairs::Vector{SVector{2,Int}}   # candidate (i<j) pairs within cutoff+skin
    r0::Vector{SVector{D,T}}        # reference positions at last rebuild
    max_disp2::T
    nbuilds::Int
    ncells::NTuple{D,Int}
    use_cells::Bool
    # reusable buffers for cell lists
    head::Vector{Int}
    next::Vector{Int}
end

function Base.show(io::IO, nl:: MyNeighborList{D,T}) where {D,T}
    print(io,
        " MyNeighborList{$D,$(T)}(",
        "pairs=", length(nl.pairs),
        ", cutoff=", nl.cutoff,
        ", skin=", nl.skin,
        ", nbuilds=", nl.nbuilds,
        ", cells=", nl.ncells,
        ", use_cells=", nl.use_cells,
        ")")
end

# --- Cell grid choice -------------------------------------------------

@inline function _choose_cell_grid(box, rskin::T) where T
    L = box.L
    return ntuple(d -> max(1, Int(floor(L / rskin))), 3)
end

# --- Construction (dispatch on ORTHO) --------------------------------

# Wrapper that receives sys::System (with sys.box)

"""
    make_neighbor_list(ps, sys; cutoff, skin, use_cells=true) ->  MyNeighborList

Build a neighbor list from the current particle positions `ps` and simulation
box `sys.box`.

Arguments
---------
- `ps::AtomicPhaseSpace{D,T}` : holds current positions `ps.r::Vector{Vec{D,T}}`
- `sys::System`               : must provide `sys.box::Box{D,T,D2,ORTHO}`
- `cutoff::T`                 : interaction cutoff (same units as positions)
- `skin::T`                   : Verlet skin buffer. Rebuild occurs when
                                `2*max_displacement_since_build > skin`.
- `use_cells::Bool=true`      : if `true` and the box type is orthorhombic
                                (`ORTHO=true`), use the O(N) linked-cell builder;
                                otherwise use the O(N²) all-pairs builder.

Details
-------
- The builder inserts i<j pairs within `rskin = cutoff + skin`. At force time,
  `iterate_pairs` filters by `r2 ≤ cutoff^2`.
- Internal buffers (`pairs`, `head`, `next`) are preallocated and reused.

Returns
-------
A ` MyNeighborList{D,T}` ready for iteration via `iterate_pairs(nl, ps, sys)`.
"""
function make_neighbor_list(ps, sys::System; cutoff::T, skin::T, use_cells::Bool=true) where {D,T<:AbstractFloat}
    _make_neighbor_list(ps, sys.box; cutoff, skin, use_cells)
end

# Main constructor specialized on Box{...,ORTHO}
function _make_neighbor_list(
    ps, box;
    cutoff::T, skin::T, use_cells::Bool=true
) where {T}
    cutoff2 = cutoff * cutoff
    rskin  = cutoff + skin
    rskin2 = rskin * rskin

    N  = length(ps.r)
    D = length(ps.r[1])
    r0 = Vector{SVector{D,T}}(undef, N)
    @inbounds for i in 1:N
        r0[i] = ps.r[i]
    end

    use_cells_here = (use_cells )
    ncells   = use_cells_here ? _choose_cell_grid(box, rskin) : ntuple(_ -> 1, D)
    ncelltot = prod(ncells)

    head = fill(0, max(ncelltot, 1))
    next = Vector{Int}(undef, N)

    pairs = Vector{SVector{2,Int}}()
    sizehint!(pairs, max(16, 100N))

    if use_cells_here
        _build_pairs_cells!(pairs, head, next, ps, box, rskin)
    else
        _build_pairs_allpairs!(pairs, ps, box, rskin)
    end

    return  MyNeighborList{D,T}(
        cutoff, cutoff2, skin, rskin2,
        pairs, r0, zero(T), 1, ncells, use_cells,
        head, next
    )
end


"""
    rebuild!(nl, ps, box)

Force a rebuild of the Verlet candidate list from the current positions.

- Reuses internal buffers (`pairs`, `head`, `next`) and only resizes if the
  particle count `N` or the cell grid shape changed.
- For orthorhombic boxes (`Box{…, ORTHO=true}`) and `nl.use_cells=true`,
  the linked-cell builder is used. For triclinic boxes (`ORTHO=false}`),
  the O(N²) all-pairs builder is used.
- Resets the displacement reference frame (`nl.r0 .= ps.r`) and `max_disp2`.

Returns `nl`.
"""
function rebuild!(nl:: MyNeighborList{D,T}, ps, box) where {D,T}
    rskin = nl.cutoff + nl.skin
    use_cells_here = (nl.use_cells)

    # Resize/clear buffers if necessary
    N = length(ps.r)
    ncells = use_cells_here ? _choose_cell_grid(box, rskin) : ntuple(_ -> 1, D)
    nct = max(prod(ncells), 1)

    if length(nl.head) != nct
        resize!(nl.head, nct)
    end
    fill!(nl.head, 0)

    if length(nl.next) != N
        resize!(nl.next, N)
    end
    empty!(nl.pairs)

    # Update grid shape
    nl.ncells = ncells

    # Rebuild candidates in-place
    if use_cells_here
        _build_pairs_cells!(nl.pairs, nl.head, nl.next, ps, box, rskin)
    else
        _build_pairs_allpairs!(nl.pairs, ps, box, rskin)
    end

    # Reset reference frame & displacement tracker
    @inbounds for i in eachindex(ps.r)
        nl.r0[i] = ps.r[i]
    end
    nl.max_disp2 = zero(T)
    nl.rskin2    = (nl.cutoff + nl.skin)^2
    nl.nbuilds  += 1
    return nl
end


fractional(box, ri) = box.L \ ri  # dispatch on ORTHO


# --- Internal builders ------------------------------------------------

# In-place cell builder — reuses buffers; uses fractional(cartesian) via dispatch.
function _build_pairs_cells!(
    pairs::Vector{SVector{2,Int}},
    head::Vector{Int},
    next::Vector{Int},
    ps,
    box,
    rskin::T
) where {T<:AbstractFloat,D2,ORTHO}
    rskin2 = rskin * rskin
    N = length(ps.r)
    D = length(ps.r[1])
    nc = _choose_cell_grid(box, rskin)
    @assert prod(nc) == length(head) "head buffer size mismatch."

    periodic = true
    fill!(head, 0)

    # Assign particles to cells (0-based per-dim indices) using fractional coords
    @inbounds for i in 1:N
        ri = ps.r[i]
        sS = fractional(box, ri)              # SVector
        s  = MVector{D,T}(undef)
        @inbounds for d in 1:D; s[d] = sS[d]; end
        _wrap_frac!(s, periodic)
        idx = ntuple(d -> begin
            sd = s[d] * nc[d]
            k  = Int(floor(sd))
            periodic ? _wrap_cell_index(k, nc[d]) : _clamp_cell_index(k, nc[d])
        end, D)
        lid = _linid0(idx, nc)
        next[i] = head[lid]
        head[lid] = i
    end

    empty!(pairs)
    # sizehint!(pairs, max(length(pairs), 100N))


    if D == 3
        @inbounds for i3 in 0:(nc[3]-1), i2 in 0:(nc[2]-1), i1 in 0:(nc[1]-1)
            lidA = _linid0((i1,i2,i3), nc)
            uniq = Vector{Int}(undef, 27); nuniq = 0
            for d3 in -1:1, d2 in -1:1, d1 in -1:1
                j1 = periodic ? _wrap_cell_index(i1 + d1, nc[1]) : (i1 + d1)
                j2 = periodic ? _wrap_cell_index(i2 + d2, nc[2]) : (i2 + d2)
                j3 = periodic ? _wrap_cell_index(i3 + d3, nc[3]) : (i3 + d3)
                (0 <= j1 < nc[1] && 0 <= j2 < nc[2] && 0 <= j3 < nc[3]) || continue
                lidB = _linid0((j1,j2,j3), nc)
                (lidB >= lidA) || continue
                seen = false
                @inbounds for t in 1:nuniq
                    if uniq[t] == lidB; seen = true; break; end
                end
                if !seen; nuniq += 1; uniq[nuniq] = lidB; end
            end
            @inbounds for t in 1:nuniq
                _accumulate_cell_pairs!(pairs, head, next, lidA, uniq[t], ps, box, rskin2)
            end
        end
    end

    return pairs
end

# Accumulate i<j pairs from cell A vs cell B, testing against rskin2
@inline function _accumulate_cell_pairs!(pairs, head, next, lidA::Int, lidB::Int, ps, box, rskin2)
    i = head[lidA]
    if lidA == lidB
        @inbounds while i != 0
            j = next[i]; ri = ps.r[i]
            while j != 0
                r2 = distance2_minimum_image(ri, ps.r[j], box)
                if r2 <= rskin2
                    push!(pairs, ifelse(i<j, SVector(i, j), SVector(j,i)))
                end
                j = next[j]
            end
            i = next[i]
        end
    else
        @inbounds while i != 0
            j = head[lidB]; ri = ps.r[i]
            while j != 0
                r2 = distance2_minimum_image(ri, ps.r[j], box)
                if r2 <= rskin2
                    if i>j
                        push!(pairs, SVector(j,i))
                    else
                        push!(pairs, SVector(i,j))
                    end
                    # push!(pairs, ifelse(i<j, SVector(i, j), SVector(j,i)))
                end
                j = next[j]
            end
            i = next[i]
        end
    end
    return nothing
end

@inline function displacement(ri::SVector{3,T}, rj::SVector{3,T}, box) where {T}
    dx = rj[1] - ri[1]
    dy = rj[2] - ri[2]
    dz = rj[3] - ri[3]
    L = box.L
    dx_corr = round(dx/L)*L
    dy_corr =  round(dy/L)*L
    dz_corr = round(dz/L)*L
    SVector{3,T}(dx - dx_corr, dy - dy_corr, dz - dz_corr)
end



@inline function displacement(ri::SVector{D,T}, rj::SVector{D,T}, box) where {D,T,D2}
    dx = rj - ri
    L = box.L
    dx_corr = @inbounds SVector{D,T}(ntuple(d ->  round(dx[d]/L)*L , D))
    dx - dx_corr
end


@inline function distance2_minimum_image(ri::SVector{D,T}, rj::SVector{D,T}, box) where {D,T,D2, ORTHO}
    dr = displacement(ri, rj, box)
    dot(dr, dr)
end




# O(N^2) fallback (triclinic or use_cells=false), in-place append to `pairs`
function _build_pairs_allpairs!(pairs::Vector{SVector{2,Int}}, ps, box, rskin::T) where {D,T<:AbstractFloat,D2,ORTHO}
    rskin2 = rskin * rskin
    N = length(ps.r)
    empty!(pairs)
    sizehint!(pairs, max(length(pairs), Int(cld(N*(N-1), 8))))
    @inbounds for i in 1:(N - 1)
        ri = ps.r[i]
        for j in (i + 1):N
            r2 = distance2_minimum_image(ri, ps.r[j], box)  # dispatches on ORTHO
            if r2 <= rskin2
                push!(pairs, SVector(i, j))
            end
        end
    end
    return pairs
end

# =====================================================================
# TEST HELPERS
# =====================================================================
"""
    allpairs_within_cutoff(ps, box, cutoff) -> Vector{SVector{2,Int}}

Reference O(N²) set of i<j pairs within `cutoff` using the minimum-image
distance in `box`. Intended for testing and validation of the cell builder.
"""
function allpairs_within_cutoff(ps, box, cutoff::T) where {D,T,D2,ORTHO}
    cutoff2 = cutoff * cutoff
    N = length(ps.r)
    out = Vector{SVector{2,Int}}()
    sizehint!(out, max(16, Int(cld(N * (N - 1), 8))))
    @inbounds for i in 1:(N - 1)
        ri = ps.r[i]
        for j in (i + 1):N
            r2 = distance2_minimum_image(ri, ps.r[j], box)
            if r2 <= cutoff2
                push!(out, SVector(i, j))
            end
        end
    end
    return out
end

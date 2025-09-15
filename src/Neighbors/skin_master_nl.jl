"""
    MasterNeighborVerletList

Master-level neighbor list that stores candidate i<j pairs within a Verlet
radius (cutoff + skin), and tracks particle displacements to decide when to
rebuild. Intended for cubic periodic boxes and 3D systems (Dims=3).

This implementation mirrors the high-level idea of the provided legacy code,
but integrates with the current package API and types.
"""
mutable struct MasterNeighborVerletList{T<:Number, IT<:Integer} <: AbstractNeighborList
    skin::T
    pairs::Vector{NTuple{2,IT}}        # candidate i<j pairs within (r_verlet + skin)
    r0::Vector{SVector{3,T}}           # reference positions at last rebuild
    max_disp2::T                       # max squared displacement since last rebuild
    grid::Union{Nothing,CellGrid{IT,T}}  # reusable cell grid
end

function MasterNeighborVerletList(skin::T; N::Integer=0) where {T}
    N = max(0, Int(N))
    pairs = NTuple{2,Int}[]
    r0    = [SVector{3,T}(0,0,0) for _ in 1:N]
    return MasterNeighborVerletList{T, Int}(skin, pairs, r0, zero(T), nothing)
end

@inline function _squared_distance_min_image(R::AbstractVector, i::Integer, j::Integer, box::CubicBox)
    Δ = R[i] - R[j]
    Δ = minimum_image(Δ, box)
    return dot(Δ, Δ)
end

# Ensure the Verlet list is valid for the given cutoff by rebuilding if needed
function _ensure_valid!(nl::MasterNeighborVerletList{T}, R::AbstractVector, box::CubicBox; r_verlet::T) where {T}
    N = length(R)
    if length(nl.r0) != N
        resize!(nl.r0, N)
        for i in 1:N
            nl.r0[i] = SVector{3,T}(R[i])
        end
        _rebuild!(nl, R, box; r_verlet)
        return nl
    end

    md2 = nl.max_disp2
    @inbounds for i in 1:N
        Δ = R[i] - nl.r0[i]
        Δ = minimum_image(Δ, box)
        d2 = dot(Δ, Δ)
        if d2 > md2
            md2 = d2
        end
    end
    nl.max_disp2 = md2
    if 2*sqrt(md2) > nl.skin
        _rebuild!(nl, R, box; r_verlet)
    end
    return nl
end

@inline function _linear_index(cx::Int, cy::Int, cz::Int, dims::NTuple{3,Int})
    nx, ny, nz = dims
    return ((cz - 1) * ny + (cy - 1)) * nx + cx
end

# Build candidate pairs within (r_verlet + skin) using a cell grid
function _rebuild!(nl::MasterNeighborVerletList{T,IT}, R::AbstractVector, box::CubicBox; r_verlet::T) where {T,IT}
    N = length(R)
    rlist = r_verlet + nl.skin
    rlist2 = rlist * rlist

    # Build or reuse cell grid sized for rlist
    g = if nl.grid === nothing || nl.grid.cell_size < rlist || nl.grid.L != box_length(box)
        build_cellgrid(R, box; cell_size=rlist)
    else
        rebin!(nl.grid, R, box)
    end
    nl.grid = g

    empty!(nl.pairs)

    nx, ny, nz = Int.(g.dims)
    heads = g.heads; nextp = g.next

    @inbounds for cz in 1:nz, cy in 1:ny, cx in 1:nx
        c_idx = _linear_index(cx, cy, cz, (nx, ny, nz))

        # Precompute unique neighbor cells for half-stencil
        neighs = StaticArrays.MVector{13,Int}(undef)
        ncount = 0
        for dz in -1:1, dy in -1:1, dx in 0:1
            if (dx == 0 && dy < 0) || (dx == 0 && dy == 0 && dz <= 0)
                continue
            end
            # Avoid double counting for very small grids
            if (nx == 2 && dx == 1 && cx != 1) ||
               (ny == 2 && dy == 1 && cy != 1) ||
               (nz == 2 && dz == 1 && cz != 1)
                continue
            end
            ncx = mod1(cx + dx, nx); ncy = mod1(cy + dy, ny); ncz = mod1(cz + dz, nz)
            cc_idx = _linear_index(ncx, ncy, ncz, (nx, ny, nz))
            if cc_idx == c_idx; continue; end
            seen = false
            @inbounds for t in 1:ncount
                if neighs[t] == cc_idx
                    seen = true; break
                end
            end
            if !seen
                ncount += 1
                neighs[ncount] = cc_idx
            end
        end

        i = heads[c_idx]
        while i != 0
            # Same-cell forward pairs
            j = nextp[i]
            while j != 0
                r2 = _squared_distance_min_image(R, i, j, box)
                if r2 <= rlist2
                    u = i < j ? i : j
                    v = i < j ? j : i
                    push!(nl.pairs, (IT(u), IT(v)))
                end
                j = nextp[j]
            end

            # Neighbor cells
            @inbounds for t in 1:ncount
                cc_idx = neighs[t]
                j = heads[cc_idx]
                while j != 0
                    r2 = _squared_distance_min_image(R, i, j, box)
                    if r2 <= rlist2
                        u = i < j ? i : j
                        v = i < j ? j : i
                        push!(nl.pairs, (IT(u), IT(v)))
                    end
                    j = nextp[j]
                end
            end

            i = nextp[i]
        end
    end

    # Reset displacement frame
    if length(nl.r0) != N
        resize!(nl.r0, N)
    end
    @inbounds for i in 1:N
        nl.r0[i] = SVector{3,T}(R[i])
    end
    nl.max_disp2 = zero(T)
    return nl
end

# --- Integration with forcefield API ----------------------------------------

function build_neighbors_from_master!(pot::AbstractPairPotential, sys::System, masterNL::MasterNeighborVerletList)
    PairType = eltype(pot.params.table)

    neighbors = pot.neighborlist.neighbors
    empty!(neighbors)

    types = sys.types
    skin  = pot.skin
    R = sys.positions
    box = sys.box

    @inbounds for (ui, uj) in masterNL.pairs
        i = Int(ui); j = Int(uj)
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

function build_all_neighbors!(master_nl::MasterNeighborVerletList, ff::ForceField, sys::System; method::Symbol=:cells)
    # Compute maximum per-layer Verlet radius (rc + skin)
    layers = ff.layers
    rc_max = maximum(
        maximum([p.rc for p in pot.params.table]) + pot.skin
        for pot in layers
    )

    if method != :cells
        error("MasterNeighborVerletList only supports method=:cells in this implementation")
    end

    # _ensure_valid!(master_nl, sys.positions, sys.box; r_verlet=rc_max)
    _rebuild!(master_nl, sys.positions, sys.box; r_verlet=rc_max)

    # Populate per-potential neighbor lists from candidate pairs
    map(pot -> build_neighbors_from_master!(pot, sys, master_nl), layers)
end


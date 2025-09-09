# Cell-based neighbor list builder (emits a half list) + lightweight support
# for using it in LJ forces & rebuild policy without altering existing APIs.
#
# Public API:
#   build_neighborlist_cells(R, box; cutoff, skin=0.3, grid=nothing) -> NeighborList-like
#
export build_neighborlist_cells

using LinearAlgebra, StaticArrays
# only pull types/builders; do not import parent API functions

# A tiny internal NL type that quacks like the existing `NeighborList` for the fields
# the tests need (.offsets and .pairs). We also carry ref_positions / cutoff / skin
# so we can implement `maybe_rebuild!` & `max_displacement_since_build` for it.
struct HalfNeighborList{T_int,T_float, Dims}
    offsets::Vector{T_int}          # CSR row offsets (length N+1)
    pairs::Vector{T_int}            # neighbor j-indices (length = number of pairs)
    ref_positions::Vector{SVector{Dims, T_float}} 
    cutoff::T_float
    skin::T_float
end



@inline function _squared_distance_min_image(R::AbstractVector, i::T_int, j::T_int, L::T_float) where {T_int,T_float}
    dr = R[i] - R[j]
    dx, dy, dz = minimum_image(dr, L)
    return dx * dx + dy * dy + dz * dz
end

@inline function _neighbors_of_cell(cx::Int, cy::Int, cz::Int, dims::NTuple{3,Int})
    nx, ny, nz = dims
    # Collect 27 neighbor cell linear indices (with periodic wrapping)
    neigh = Vector{Int}(undef, 27)
    k = 1
    @inbounds for dz in (-1:1), dy in (-1:1), dx in (-1:1)
        ncx = mod1(cx + dx, nx)   # 1..nx
        ncy = mod1(cy + dy, ny)   # 1..ny
        ncz = mod1(cz + dz, nz)   # 1..nz
        neigh[k] = ((ncz - 1) * ny + (ncy - 1)) * nx + ncx
        k += 1
    end
    return neigh
end

"""
    build_neighborlist_cells(R, box; cutoff, skin=0.3, grid=nothing) -> NeighborList-like

Construct a **half neighbor list** (each pair appears exactly once with j>i) using
an O(N) cell-linked grid build. Returns an object with `.offsets` and `.pairs`
fields compatible with typical CSR traversal.
"""
function build_neighborlist_cells(R::AbstractVector, box;
    cutoff::Real,
    skin::Real=0.3,
    grid::Union{Nothing,CellGrid}=nothing)
    @assert length(R[1]) == 3 "d=3 only"
    N = length(R)
    L = float(box_length(box))
    rlist = float(cutoff + skin)
    rlist2 = rlist * rlist

    # Ensure we have a grid with cell_size ≥ rlist (effective cell width must be at least rlist)
    g = if grid === nothing || grid.cell_size < rlist || grid.L != L
        build_cellgrid(R, box; cell_size=rlist)
    else
        rebin!(grid, R, box)
    end

    nx, ny, nz = Int.(g.dims)
    # Precompute head pointers for convenience
    heads = g.heads
    nextp = g.next

    # First pass: count neighbors per "owner" row (owner = min(i, j))
    count_per_i = zeros(Int, N)
    @inbounds for cz in 1:nz, cy in 1:ny, cx in 1:nx
        c = ((cz - 1) * ny + (cy - 1)) * nx + cx
        # Precompute neighbor cell ids once
        neigh = _neighbors_of_cell(cx, cy, cz, (nx, ny, nz))
        # iterate particles i in cell c
        i = heads[c]
        while i != 0
            # Same-cell: only j after i (half, j>i)
            j = nextp[i]
            while j != 0
                r2 = _squared_distance_min_image(R, i, j, L)
                if r2 <= rlist2 + 1e-12
                    owner = ifelse(i < j, i, j)  # ensure row owner=min(i,j)
                    count_per_i[owner] += 1
                end
                j = nextp[j]
            end
            # Cross-cell: only consider neighbor cells with cc > c
            for idx in 1:27
                cc = neigh[idx]
                if cc <= c
                    continue
                end
                j = heads[cc]
                while j != 0
                    r2 = _squared_distance_min_image(R, i, j, L)
                    if r2 <= rlist2 + 1e-12
                        owner = ifelse(i < j, i, j)
                        count_per_i[owner] += 1
                    end
                    j = nextp[j]
                end
            end
            i = nextp[i]
        end
    end

    # Build CSR offsets
    offsets = Vector{Int}(undef, N + 1)
    offsets[1] = 1
    @inbounds for i in 1:N
        offsets[i+1] = offsets[i] + count_per_i[i]
    end
    total_pairs = offsets[end] - 1
    pairs = Vector{Int}(undef, total_pairs)

    # Second pass: fill pairs
    write_ptr = copy(offsets)
    @inbounds for cz in 1:nz, cy in 1:ny, cx in 1:nx
        c = ((cz - 1) * ny + (cy - 1)) * nx + cx
        neigh = _neighbors_of_cell(cx, cy, cz, (nx, ny, nz))
        i = heads[c]
        while i != 0
            # Same-cell: j after i only
            j = nextp[i]
            while j != 0
                r2 = _squared_distance_min_image(R, i, j, L)
                if r2 <= rlist2 + 1e-12
                    owner = ifelse(i < j, i, j)
                    other = ifelse(i < j, j, i)   # ensure other > owner
                    w = write_ptr[owner]
                    pairs[w] = other
                    write_ptr[owner] = w + 1
                end
                j = nextp[j]
            end
            # Cross-cell: only cc > c
            for idx in 1:27
                cc = neigh[idx]
                if cc <= c
                    continue
                end
                j = heads[cc]
                while j != 0
                    r2 = _squared_distance_min_image(R, i, j, L)
                    if r2 <= rlist2 + 1e-12
                        owner = ifelse(i < j, i, j)
                        other = ifelse(i < j, j, i)
                        w = write_ptr[owner]
                        pairs[w] = other
                        write_ptr[owner] = w + 1
                    end
                    j = nextp[j]
                end
            end
            i = nextp[i]
        end
    end

    return HalfNeighborList(offsets, pairs, copy(R), float(cutoff), float(skin))
end


# Use Verlet.minimum_image! to mirror the standard NeighborList logic exactly
function max_displacement_since_build(nl::Union{HalfNeighborList, NeighborList},
    R::AbstractVector, box)
    N = length(R)
    D = length(R[1])
    # Docstring is attached to the public NeighborList method; here we mirror
    # the semantics so that `maybe_rebuild!` works identically for cell lists.
    # We use the same minimum-image convention as in the rest of the package.
    @assert size(nl.ref_positions) == (N,)
    @assert length(nl.ref_positions[1]) == D
    Δ = zeros(eltype(nl.ref_positions))
    maxd2 = 0.0
    @inbounds for i in 1:N
        Δ = minimum_image(R[i] - nl.ref_positions[i], box)
        # critical: same minimum-image convention as the rest of the package
        d2 = dot(Δ, Δ)
        if d2 > maxd2
            maxd2 = d2
        end
    end
    return sqrt(maxd2)
end

function maybe_rebuild!(nl::HalfNeighborList,
    R::AbstractVector, box)
    # Classic half-skin rule (rebuild when max displacement > skin/2).
    # With O(N) builds you can safely choose a smaller skin (e.g. 0.2–0.3)
    # to tighten force errors without paying quadratic rebuild costs.
    return max_displacement_since_build(nl, R, box) > (nl.skin / 2)
end

# LJ using CSR half neighbor list (j>i)
function lj_forces(R::AbstractVector, box, nl::HalfNeighborList;
    ϵ::Real=1.0, σ::Real=1.0, rcut::Real=2.5,
    shift::Bool=true, return_potential::Bool=false)
    # Half-list aware kernel: each pair is visited once (j>i guaranteed).
    # This removes a branch and halves memory traffic relative to a full list.
    # Potential can be shifted so that U(rcut)=0 (no force smoothing).
    N = size(R, 1)
    L = T_float(box_length(box))
    rcut2 = T_float(rcut)^2
    σ2 = T_float(σ)^2
    F = zeros(eltype(R), N)

    Ushift = 0.0
    if shift
        invr2c = 1.0 / rcut2
        s2c = σ2 * invr2c
        s6c = s2c^3
        Ushift = 4.0 * float(ϵ) * (s6c^2 - s6c)
    end

    Uacc = 0.0
    @inbounds for i in 1:N
        for idx in nl.offsets[i]:(nl.offsets[i+1]-1)
            j = nl.pairs[idx]   # j>i
            dr = R[i] - R[j]
            Δ = minimum_image(dr, L)
            r2 = dot(Δ, Δ)
            if r2 <= rcut2
                invr2 = 1.0 / r2
                s2 = σ2 * invr2
                s6 = s2^3
                fr_over_r = 24.0 * float(ϵ) * (2.0 * s6^2 - s6) * invr2
                f = fr_over_r * Δ
                
                F[i] += f
                F[j] -= f 
                if return_potential
                    Uacc += 4.0 * float(ϵ) * (s6^2 - s6) - (shift ? Ushift : 0.0)
                end
            end
        end
    end
    return return_potential ? (F, Uacc) : F
end


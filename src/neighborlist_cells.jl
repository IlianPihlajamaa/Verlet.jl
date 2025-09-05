# Cell-based neighbor list builder (emits a half list) + lightweight support
# for using it in LJ forces & rebuild policy without altering existing APIs.
#
# Public API:
#   build_neighborlist_cells(R, box; cutoff, skin=0.3, grid=nothing) -> NeighborList-like
#
module __NLCellsInternal__
export build_neighborlist_cells

using LinearAlgebra
# only pull types/builders; do not import parent API functions
using ..Verlet: CellGrid, build_cellgrid, rebin!

# A tiny internal NL type that quacks like the existing `NeighborList` for the fields
# the tests need (.offsets and .pairs). We also carry ref_positions / cutoff / skin
# so we can implement `maybe_rebuild!` & `max_displacement_since_build` for it.
struct HalfNeighborList
    offsets::Vector{Int}          # CSR row offsets (length N+1)
    pairs::Vector{Int}            # neighbor j-indices (length = number of pairs)
    ref_positions::Matrix{Float64}
    cutoff::Float64
    skin::Float64
end

@inline function _box_length(box)
    try
        return getfield(box, :L)
    catch
        return getproperty(box, :L)
    end
end

@inline function _min_image!(Δ::NTuple{3,Float64}, L::Float64)
    dx, dy, dz = Δ
    dx -= L * round(dx / L)
    dy -= L * round(dy / L)
    dz -= L * round(dz / L)
    return (dx, dy, dz)
end

@inline function _squared_distance_min_image(R::AbstractMatrix, i::Int, j::Int, L::Float64)
    dx = Float64(R[i,1]) - Float64(R[j,1])
    dy = Float64(R[i,2]) - Float64(R[j,2])
    dz = Float64(R[i,3]) - Float64(R[j,3])
    dx, dy, dz = _min_image!((dx,dy,dz), L)
    return dx*dx + dy*dy + dz*dz
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
function build_neighborlist_cells(R::AbstractMatrix, box;
                                  cutoff::Real,
                                  skin::Real=0.3,
                                  grid::Union{Nothing,CellGrid}=nothing)
    @assert size(R,2) ≥ 3 "R must be N×3 (at least 3 columns)."
    N = size(R,1)
    L = float(_box_length(box))
    rlist = float(cutoff + skin)
    rlist2 = rlist*rlist

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
    offsets = Vector{Int}(undef, N+1)
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

    return HalfNeighborList(offsets, pairs, copy(Array{Float64}(R)), float(cutoff), float(skin))
end

# === Attach HalfNeighborList methods to the parent module (Verlet) ==========

const _PARENT = Base.parentmodule(@__MODULE__)

@eval _PARENT begin
    # Use Verlet.minimum_image! to mirror the standard NeighborList logic exactly
    function max_displacement_since_build(nl::__NLCellsInternal__.HalfNeighborList,
                                          R::AbstractMatrix, box)
        N, D = size(R)
    # Docstring is attached to the public NeighborList method; here we mirror
    # the semantics so that `maybe_rebuild!` works identically for cell lists.
    # We use the same minimum-image convention as in the rest of the package.
        @assert size(nl.ref_positions) == (N, D)
        Δ = zeros(eltype(nl.ref_positions), D)
        maxd2 = 0.0
        @inbounds for i in 1:N
            @inbounds for k in 1:D
                Δ[k] = R[i,k] - nl.ref_positions[i,k]
            end
            # critical: same minimum-image convention as the rest of the package
            minimum_image!(Δ, box)
            d2 = 0.0
            @inbounds for k in 1:D
                d2 += Δ[k]*Δ[k]
            end
            if d2 > maxd2
                maxd2 = d2
            end
        end
        return sqrt(maxd2)
    end

    function maybe_rebuild!(nl::__NLCellsInternal__.HalfNeighborList,
                            R::AbstractMatrix, box)
    # Classic half-skin rule (rebuild when max displacement > skin/2).
    # With O(N) builds you can safely choose a smaller skin (e.g. 0.2–0.3)
    # to tighten force errors without paying quadratic rebuild costs.
        return max_displacement_since_build(nl, R, box) > (nl.skin / 2)
    end

    # LJ using CSR half neighbor list (j>i)
    function lj_forces(R::AbstractMatrix, box, nl::__NLCellsInternal__.HalfNeighborList;
                       ϵ::Real=1.0, σ::Real=1.0, rcut::Real=2.5,
                       shift::Bool=true, return_potential::Bool=false)
    # Half-list aware kernel: each pair is visited once (j>i guaranteed).
    # This removes a branch and halves memory traffic relative to a full list.
    # Potential can be shifted so that U(rcut)=0 (no force smoothing).
        N = size(R,1)
        L = float(__NLCellsInternal__._box_length(box))
        rcut2 = float(rcut)^2
        σ2 = float(σ)^2
        F = zeros(Float64, N, 3)

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
                dx = Float64(R[i,1]) - Float64(R[j,1])
                dy = Float64(R[i,2]) - Float64(R[j,2])
                dz = Float64(R[i,3]) - Float64(R[j,3])
                dx -= L * round(dx / L)
                dy -= L * round(dy / L)
                dz -= L * round(dz / L)
                r2 = dx*dx + dy*dy + dz*dz
                if r2 == 0.0
                    continue
                end
                if r2 <= rcut2
                    invr2 = 1.0 / r2
                    s2 = σ2 * invr2
                    s6 = s2^3
                    fr_over_r = 24.0 * float(ϵ) * (2.0*s6^2 - s6) * invr2
                    fx = fr_over_r * dx
                    fy = fr_over_r * dy
                    fz = fr_over_r * dz
                    F[i,1] += fx; F[i,2] += fy; F[i,3] += fz
                    F[j,1] -= fx; F[j,2] -= fy; F[j,3] -= fz
                    if return_potential
                        Uacc += 4.0 * float(ϵ) * (s6^2 - s6) - (shift ? Ushift : 0.0)
                    end
                end
            end
        end
        return return_potential ? (F, Uacc) : F
    end
end

end # module __NLCellsInternal__

using .__NLCellsInternal__: build_neighborlist_cells
export build_neighborlist_cells

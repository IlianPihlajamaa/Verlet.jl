
"""
    CellGrid{D,IT<:Integer,T<:Real}

Linked-list **cell grid** for cubic periodic boxes. The domain is split into a
uniform `dims` grid of cubic cells with an **intrusive linked list** (`heads`,
`next`) storing particle indices in each cell. This enables **O(N)** rebinning
at fixed density and powers the O(N) neighbor build.

# Fields
- `L::T`: cubic box length.
- `cell_size::T`: effective uniform cell width used for binning.
- `dims::NTuple{D,IT}`: number of cells along each axis (each ≥ 1).
- `heads::Vector{IT}`: length `prod(dims)`; head index per cell (`0` sentinel).
- `next::Vector{IT}`: length `N`; linked list “next” pointer per particle (`0`).

!!! tip "How `cell_size` is chosen"
    `build_cellgrid` computes `n = floor(Int, L/cell_size)` and then uses the
    **effective** width `L/n` for indexing so that a `(3ᵈ)` stencil is sufficient
    for a search radius ≤ `cell_size`.

!!! warning "Units"
    `R` and `L` must be expressed in the **same length units**.
"""
struct CellGrid{D,IT,T}
    L::T
    cell_size::T
    dims::NTuple{D,IT}
    heads::Vector{IT}
    next::Vector{IT}
end



@inline _n_cells(L::Real, cell_size::Real) = max(1, Int(floor(L / cell_size)))

@inline function _linear_index(indices::NTuple{D,Int}, dims::NTuple{D,Int}) where {D}
    id = 1
    stride = 1
    @inbounds for d in 1:D
        id += (indices[d] - 1) * stride
        stride *= dims[d]
    end
    return id
end

@inline function _coord_to_cell_idx(x::T_Float, L::T_Float, cell_size::T_Float, n::Int)
    # Map position to [0, L) robustly, then to 1..n
    x0 = x + 0.5 * L
    x0 -= floor(x0 / L) * L   # now in [0, L)
    c = Int(floor(x0 / cell_size)) + 1
    return c > n ? n : (c < 1 ? 1 : c)
end

# -- Build & Rebin -------------------------------------------------------------

"""
    build_cellgrid(R, box; cell_size)

Create a new `CellGrid` sized for `cell_size` and bin positions `R` (N×D).
Returns a populated grid with linked lists set for particle indices 1..N.
"""
function build_cellgrid(R::AbstractVector, box; cell_size::Real)
    D = length(R[1])
    return _build_cellgrid(R, box, cell_size, Val(D))
end

function _build_cellgrid(R::AbstractVector, box, cell_size::Real, ::Val{D}) where {D}
    @assert length(R[1]) == D "d=$D"
    N = length(R)
    L = float(box_length(box))
    n = _n_cells(L, float(cell_size))
    cs_eff = L / n
    dims = ntuple(_ -> n, D)
    heads = fill(Int(0), max(prod(dims), 1))
    nxt = fill(Int(0), N)
    grid = CellGrid{D,Int,T_Float}(L, cs_eff, dims, heads, nxt)
    return rebin!(grid, R, box)
end

"""
    rebin!(grid, R, box) -> grid

Reset the grid's `heads` and `next` and bin the positions `R` into cells
according to the current `cell_size` and periodic cubic box `box`.
"""
function rebin!(grid::CellGrid{D,IT,T}, R::AbstractVector, box) where {D,IT<:Integer,T<:Real}
    @assert length(R[1]) == D "d=$D"
    N = length(R)
    L = T_Float(grid.L)
    dims = ntuple(d -> Int(grid.dims[d]), D)
    cs = T_Float(grid.cell_size)

    fill!(grid.heads, IT(0))
    if length(grid.next) != N
        resize!(grid.next, N)
    end
    fill!(grid.next, IT(0))

    @inbounds for i in 1:N
        r = R[i]
        idxs = ntuple(d -> _coord_to_cell_idx(T_Float(r[d]), L, cs, dims[d]), D)
        c = _linear_index(idxs, dims)
        grid.next[i] = grid.heads[c]
        grid.heads[c] = IT(i)
    end
    return grid
end

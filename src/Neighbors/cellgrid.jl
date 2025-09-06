
"""
    CellGrid{IT<:Integer, T<:Real}

Linked-list **cell grid** for cubic periodic boxes. The domain is split into a
uniform `dims = (nx,ny,nz)` grid of cubic cells with an **intrusive linked list**
(`heads`, `next`) storing particle indices in each cell. This enables **O(N)**
rebinning at fixed density and powers the O(N) neighbor build.

# Fields
- `L::T`: cubic box length.
- `cell_size::T`: **effective** uniform cell width actually used for binning.
- `dims::NTuple{3,IT}`: number of cells along each axis (each ≥ 1).
- `heads::Vector{IT}`: length `nx*ny*nz`; head index per cell (`0` sentinel = empty).
- `next::Vector{IT}`: length `N`; linked list “next” pointer per particle (`0` = end).

!!! tip "How `cell_size` is chosen"
    `build_cellgrid` computes `nx = floor(Int, L/cell_size)` and then uses the
    **effective** width `L/nx` for indexing so that a 27-cell sweep is sufficient
    for a search radius ≤ `cell_size`.

!!! warning "Units"
    `R` and `L` must be expressed in the **same length units**.
"""
struct CellGrid{IT,T}
    L::T                  # cubic box length
    cell_size::T          # typical choice: cutoff + skin (rlist)
    dims::NTuple{3,IT}    # (nx, ny, nz), each ≥ 1
    heads::Vector{IT}     # length = nx*ny*nz; 0 sentinel = empty
    next::Vector{IT}      # length = N; next particle index in the cell list (0 = end)
end



@inline _n_cells(L::Real, cell_size::Real) = max(1, Int(floor(L / cell_size)))

@inline function _linear_index(cx::Int, cy::Int, cz::Int, dims::NTuple{3,Int})
    nx, ny, nz = dims
    # 1-based linearization
    return ((cz - 1) * ny + (cy - 1)) * nx + cx
end

@inline function _coord_to_cell_idx(x::Float64, L::Float64, cell_size::Float64, n::Int)
    # Map position to [0, L) robustly, then to 1..n
    x0 = x + 0.5 * L
    x0 -= floor(x0 / L) * L   # now in [0, L)
    c = Int(floor(x0 / cell_size)) + 1
    return c > n ? n : (c < 1 ? 1 : c)
end

# -- Build & Rebin -------------------------------------------------------------

"""
    build_cellgrid(R, box; cell_size)

Create a new `CellGrid` sized for `cell_size` and bin positions `R` (N×3).
Returns a populated grid with linked lists set for particle indices 1..N.
"""
function build_cellgrid(R::AbstractMatrix, box; cell_size::Real)
    @assert size(R,2) ≥ 3 "R must be N×3 (at least 3 columns)."
    N = size(R,1)
    L = float(box_length(box))
    # Choose number of cells so that the EFFECTIVE uniform width cs_eff = L/n ≥ requested cell_size.
    # This guarantees the standard 27-neighbor sweep is sufficient for a search radius ≤ cell_size.
    nx = _n_cells(L, float(cell_size))              # floor(L / requested)
    cs_eff = L / nx                                 # effective width (≥ requested)
    dims = (nx, nx, nx)
    heads = fill(Int(0), nx*nx*nx)
    nxt = fill(Int(0), N)
    # Store the *effective* cell size actually used for binning (uniform tiling)
    grid = CellGrid{Int,Float64}(L, cs_eff, dims, heads, nxt)
    return rebin!(grid, R, box)
end

"""
    rebin!(grid, R, box) -> grid

Reset the grid's `heads` and `next` and bin the positions `R` into cells
according to the current `cell_size` and periodic cubic box `box`.
"""
function rebin!(grid::CellGrid{IT,T}, R::AbstractMatrix, box) where {IT<:Integer,T<:Real}
    @assert size(R,2) ≥ 3 "R must be N×3 (at least 3 columns)."
    N = size(R,1)
    L = Float64(grid.L)
    nx, ny, nz = Int.(grid.dims)
    cs = Float64(grid.cell_size)  # effective uniform width used for indexing

    # Reset lists
    fill!(grid.heads, IT(0))
    if length(grid.next) != N
        resize!(grid.next, N)
    end
    fill!(grid.next, IT(0))

    @inbounds for i in 1:N
        cx = _coord_to_cell_idx(Float64(R[i,1]), L, cs, nx)
        cy = _coord_to_cell_idx(Float64(R[i,2]), L, cs, ny)
        cz = _coord_to_cell_idx(Float64(R[i,3]), L, cs, nz)
        c  = _linear_index(cx, cy, cz, (nx, ny, nz))
        # push-front i into cell c
        grid.next[i] = grid.heads[c]
        grid.heads[c] = IT(i)
    end
    return grid
end

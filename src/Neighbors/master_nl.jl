<<<<<<< HEAD
using LinearAlgebra, StaticArrays

@inline function _squared_distance_min_image(R::AbstractVector, i::T_int, j::T_int, box::CubicBox)
=======
@inline function _squared_distance_min_image(R::AbstractVector, i::T_Int, j::T_Int, box::CubicBox)
>>>>>>> e020923 (fix all namespaces)
    dr = R[i] - R[j]
    dr = minimum_image(dr, box)
    return dot(dr, dr)
end

function _build_master_nl_bruteforce!(entries, R, box, rlist2)
    N = length(R)
    @inbounds for i in 1:N-1
        for j in i+1:N
            r2 = _squared_distance_min_image(R, i, j, box)
            if r2 <= rlist2
                push!(entries, MasterNeighborEntry(i, j, r2))
            end
        end
    end
end

function _build_master_nl_allpairs!(entries, R, box)
    N = length(R)
    @inbounds for i in 1:N-1
        for j in i+1:N
            r2 = _squared_distance_min_image(R, i, j, box)
            push!(entries, MasterNeighborEntry(i, j, r2))
        end
    end
end

@inline function _neighbors_of_cell(cx::Int, cy::Int, cz::Int, dims::NTuple{3,Int})
    nx, ny, nz = dims
    neigh = Vector{Int}(undef, 27)
    k = 1
    @inbounds for dz in (-1:1), dy in (-1:1), dx in (-1:1)
        ncx = mod1(cx + dx, nx)
        ncy = mod1(cy + dy, ny)
        ncz = mod1(cz + dz, nz)
        neigh[k] = ((ncz - 1) * ny + (ncy - 1)) * nx + ncx
        k += 1
    end
    return neigh
end

function _build_master_nl_cells!(entries, R, box, rlist2, grid)
    N = length(R)
    rlist = sqrt(rlist2)
    g = if grid === nothing || grid.cell_size < rlist || grid.L != box.L
        build_cellgrid(R, box; cell_size=rlist)
    else
        rebin!(grid, R, box)
    end

    nx, ny, nz = Int.(g.dims)
    heads = g.heads
    nextp = g.next

    @inbounds for cz in 1:nz, cy in 1:ny, cx in 1:nx
        c_idx = ((cz - 1) * ny + (cy - 1)) * nx + cx
        i = heads[c_idx]
        while i != 0
            # Same cell pairs
            j = nextp[i]
            while j != 0
                r2 = _squared_distance_min_image(R, i, j, box)
                if r2 <= rlist2
                    push!(entries, MasterNeighborEntry(min(i,j), max(i,j), r2))
                end
                j = nextp[j]
            end

            # Neighbor cell pairs (27 neighbors)
            neigh = _neighbors_of_cell(cx, cy, cz, (nx,ny,nz))
            for cc_idx in neigh
                if c_idx < cc_idx
                    j = heads[cc_idx]
                    while j != 0
                        r2 = _squared_distance_min_image(R, i, j, box)
                        if r2 <= rlist2
                            push!(entries, MasterNeighborEntry(min(i,j), max(i,j), r2))
                        end
                        j = nextp[j]
                    end
                end
            end
            i = nextp[i]
        end
    end
end


"""
    build_master_neighborlist!(master_nl, R, box; r_verlet, skin, method=:cells, grid=nothing)

Update a master neighbor list in-place.

The `method` keyword argument can be one of:
- `:cells`: Use an O(N) cell-linked grid build (default).
- `:bruteforce`: Use an O(N^2) brute-force build.
- `:all_pairs`: Include all pairs, ignoring the cutoff.
"""
function build_master_neighborlist!(master_nl::MasterNeighborList, R::AbstractVector, box::CubicBox;
    r_verlet::Real=0.0,
    method::Symbol=:cells,
    grid::Union{Nothing,CellGrid}=nothing)

    @assert length(R[1]) == 3 "d=3 only"
    resize!(master_nl.entries, 0)

    if method == :cells
        rlist2 = (r_verlet + master_nl.skin)^2
        _build_master_nl_cells!(master_nl.entries, R, box, rlist2, grid)
    elseif method == :bruteforce
        rlist2 = (r_verlet + master_nl.skin)^2
        _build_master_nl_bruteforce!(master_nl.entries, R, box, rlist2)
    elseif method == :all_pairs
        _build_master_nl_allpairs!(master_nl.entries, R, box)
    else
        error("Unknown neighborlist method: $method")
    end

    if method == :cells
        sort!(master_nl.entries, by = x -> (x.i, x.j))
        unique!(master_nl.entries)
    end

    return master_nl
end

@inline function _squared_distance_min_image(R::AbstractVector, i::T_Int, j::T_Int, box::CubicBox)
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

        # Precompute unique neighbor cells for this base cell using a 13-cell half-stencil
        neighs = StaticArrays.MVector{13,Int}(undef)
        ncount = 0
        for dz in -1:1
            for dy in -1:1
                for dx in 0:1
                    # skip non-forward directions and the origin
                    if dx == 0 && dy < 0
                        continue
                    end
                    if dx == 0 && dy == 0 && dz <= 0
                        continue
                    end
                    # Tie-breakers for grids with n==2 along any axis to prevent double counting
                    if (nx == 2 && dx == 1 && cx != 1) ||
                       (ny == 2 && dy == 1 && cy != 1) ||
                       (nz == 2 && dz == 1 && cz != 1)
                        continue
                    end
                    ncx = mod1(cx + dx, nx)
                    ncy = mod1(cy + dy, ny)
                    ncz = mod1(cz + dz, nz)
                    cc_idx = ((ncz - 1) * ny + (ncy - 1)) * nx + ncx
                    # Skip self and deduplicate offsets that wrap to same cell
                    if cc_idx == c_idx
                        continue
                    end
                    seen = false
                    @inbounds for t in 1:ncount
                        if neighs[t] == cc_idx
                            seen = true
                            break
                        end
                    end
                    if !seen
                        ncount += 1
                        neighs[ncount] = cc_idx
                    end
                end
            end
        end

        i = heads[c_idx]
        while i != 0
            # Same-cell pairs: traverse forward in the intrusive list; enforce i<j ordering
            j = nextp[i]
            while j != 0
                r2 = _squared_distance_min_image(R, i, j, box)
                if r2 <= rlist2
                    if i < j
                        push!(entries, MasterNeighborEntry(i, j, r2))
                    else
                        push!(entries, MasterNeighborEntry(j, i, r2))
                    end
                end
                j = nextp[j]
            end

            # Neighbor cell pairs: iterate precomputed unique neighbor cells
            @inbounds for t in 1:ncount
                cc_idx = neighs[t]
                j = heads[cc_idx]
                while j != 0
                    r2 = _squared_distance_min_image(R, i, j, box)
                    if r2 <= rlist2
                        # Enforce i<j invariant expected elsewhere
                        if i < j
                            entry = MasterNeighborEntry(i, j, r2)
                            push!(entries, entry)
                        else
                            push!(entries, MasterNeighborEntry(j, i, r2))
                        end
                    end
                    j = nextp[j]
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
    empty!(master_nl.entries)

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

    # if method == :cells
    #     sort!(master_nl.entries, by = x -> (x.i, x.j))
    #     unique!(master_nl.entries)
    # end

    return master_nl
end

 
##############
# Preallocated two-pass CSR builder helpers
##############

@inline function _count_pairs_cells!(counts::Vector{<:Integer}, R, box, rlist2, g)
    nx, ny, nz = Int.(g.dims)
    heads = g.heads; nextp = g.next
    @inbounds for cz in 1:nz, cy in 1:ny, cx in 1:nx
        c_idx = ((cz - 1) * ny + (cy - 1)) * nx + cx

        # Precompute unique neighbor cells (half-stencil)
        neighs = StaticArrays.MVector{13,Int}(undef); ncount = 0
        for dz in -1:1, dy in -1:1, dx in 0:1
            if (dx == 0 && dy < 0) || (dx == 0 && dy == 0 && dz <= 0)
                continue
            end
            ncx = mod1(cx + dx, nx); ncy = mod1(cy + dy, ny); ncz = mod1(cz + dz, nz)
            cc_idx = ((ncz - 1) * ny + (ncy - 1)) * nx + ncx
            if cc_idx == c_idx; continue; end
            seen = false
            @inbounds for t in 1:ncount
                if neighs[t] == cc_idx; seen = true; break; end
            end
            if !seen
                ncount += 1; neighs[ncount] = cc_idx
            end
        end

        i = heads[c_idx]
        while i != 0
            # Same-cell
            j = nextp[i]
            while j != 0
                r2 = _squared_distance_min_image(R, i, j, box)
                if r2 <= rlist2
                    u = i < j ? i : j
                    counts[u] += 1
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
                        counts[u] += 1
                    end
                    j = nextp[j]
                end
            end
            i = nextp[i]
        end
    end
end

@inline function _fill_csr_cells!(rowptr, writeptr, colind, R, box, rlist2, g)
    nx, ny, nz = Int.(g.dims)
    heads = g.heads; nextp = g.next
    copyto!(writeptr, rowptr)  # running write positions (1-based)
    @inbounds for cz in 1:nz, cy in 1:ny, cx in 1:nx
        c_idx = ((cz - 1) * ny + (cy - 1)) * nx + cx
        neighs = StaticArrays.MVector{13,Int}(undef); ncount = 0
        for dz in -1:1, dy in -1:1, dx in 0:1
            if (dx == 0 && dy < 0) || (dx == 0 && dy == 0 && dz <= 0)
                continue
            end
            ncx = mod1(cx + dx, nx); ncy = mod1(cy + dy, ny); ncz = mod1(cz + dz, nz)
            cc_idx = ((ncz - 1) * ny + (ncy - 1)) * nx + ncx
            if cc_idx == c_idx; continue; end
            seen = false
            @inbounds for t in 1:ncount
                if neighs[t] == cc_idx; seen = true; break; end
            end
            if !seen
                ncount += 1; neighs[ncount] = cc_idx
            end
        end

        i = heads[c_idx]
        while i != 0
            j = nextp[i]
            while j != 0
                r2 = _squared_distance_min_image(R, i, j, box)
                if r2 <= rlist2
                    if i < j
                        k = writeptr[i]; colind[k] = j; writeptr[i] = k + 1
                    else
                        k = writeptr[j]; colind[k] = i; writeptr[j] = k + 1
                    end
                end
                j = nextp[j]
            end
            @inbounds for t in 1:ncount
                cc_idx = neighs[t]
                if cc_idx <= c_idx
                    continue
                end
                j = heads[cc_idx]
                while j != 0
                    r2 = _squared_distance_min_image(R, i, j, box)
                    if r2 <= rlist2
                        if i < j
                            k = writeptr[i]; colind[k] = j; writeptr[i] = k + 1
                        else
                            k = writeptr[j]; colind[k] = i; writeptr[j] = k + 1
                        end
                    end
                    j = nextp[j]
                end
            end
            i = nextp[i]
        end
    end
end

@inline function _count_pairs_bruteforce!(counts, R, box, rlist2)
    N = length(R)
    @inbounds for i in 1:N-1
        for j in i+1:N
            r2 = _squared_distance_min_image(R, i, j, box)
            if r2 <= rlist2
                counts[i] += 1
            end
        end
    end
end

@inline function _fill_csr_bruteforce!(rowptr, writeptr, colind, R, box, rlist2)
    N = length(R)
    copyto!(writeptr, rowptr)
    @inbounds for i in 1:N-1
        for j in i+1:N
            r2 = _squared_distance_min_image(R, i, j, box)
            if r2 <= rlist2
                k = writeptr[i]; colind[k] = j; writeptr[i] = k + 1
            end
        end
    end
end

@inline function _count_pairs_allpairs!(counts, R)
    N = length(R)
    @inbounds for i in 1:N-1
        counts[i] += (N - i)
    end
end

@inline function _fill_csr_allpairs!(rowptr, writeptr, colind, R)
    N = length(R)
    copyto!(writeptr, rowptr)
    @inbounds for i in 1:N-1
        for j in i+1:N
            k = writeptr[i]; colind[k] = j; writeptr[i] = k + 1
        end
    end
end

function build_master_neighborlist!(master_nl::MasterNeighborCSRList, R::AbstractVector, box::CubicBox;
    r_verlet::Real=0.0, method::Symbol=:cells, grid::Union{Nothing,CellGrid}=nothing)

    @assert length(R[1]) == 3 "d=3 only"
    N = length(R)

    # Ensure preallocated buffers have correct sizes
    if length(master_nl.rowptr) != N + 1
        resize!(master_nl.rowptr, N + 1)
    end
    if length(master_nl.counts) != N
        resize!(master_nl.counts, N)
    end
    if length(master_nl.writeptr) != N + 1
        resize!(master_nl.writeptr, N + 1)
    end
    fill!(master_nl.counts, 0)

    # Prepare per-row scratch buffers (reused)
    rows = master_nl.rows
    if length(rows) != N
        resize!(rows, N)
        for i in 1:N
            rows[i] = Int[]
        end
    end
    @inbounds for i in 1:N
        empty!(rows[i])
    end

    # Single pass: compute distances once, append to row buffers
    if method == :cells
        rlist2 = (r_verlet + master_nl.skin)^2
        rlist = sqrt(rlist2)
        g = if grid === nothing || grid.cell_size < rlist || grid.L != box.L
            build_cellgrid(R, box; cell_size=rlist)
        else
            rebin!(grid, R, box)
        end

        nx, ny, nz = Int.(g.dims)
        heads = g.heads; nextp = g.next
        @inbounds for cz in 1:nz, cy in 1:ny, cx in 1:nx
            c_idx = ((cz - 1) * ny + (cy - 1)) * nx + cx
            neighs = StaticArrays.MVector{13,Int}(undef); ncount = 0
            for dz in -1:1, dy in -1:1, dx in 0:1
                if (dx == 0 && dy < 0) || (dx == 0 && dy == 0 && dz <= 0)
                    continue
                end
                # Tie-breakers for n==2 along any axis
                if (nx == 2 && dx == 1 && cx != 1) ||
                   (ny == 2 && dy == 1 && cy != 1) ||
                   (nz == 2 && dz == 1 && cz != 1)
                    continue
                end
                ncx = mod1(cx + dx, nx); ncy = mod1(cy + dy, ny); ncz = mod1(cz + dz, nz)
                cc_idx = ((ncz - 1) * ny + (ncy - 1)) * nx + ncx
                if cc_idx == c_idx; continue; end
                seen = false
                @inbounds for t in 1:ncount
                    if neighs[t] == cc_idx; seen = true; break; end
                end
                if !seen
                    ncount += 1; neighs[ncount] = cc_idx
                end
            end

            i = heads[c_idx]
            while i != 0
                j = nextp[i]
                while j != 0
                    r2 = _squared_distance_min_image(R, i, j, box)
                    if r2 <= rlist2
                        u = i < j ? i : j
                        v = i < j ? j : i
                        push!(rows[u], v)
                    end
                    j = nextp[j]
                end
                @inbounds for t in 1:ncount
                    cc_idx = neighs[t]
                    j = heads[cc_idx]
                    while j != 0
                        r2 = _squared_distance_min_image(R, i, j, box)
                        if r2 <= rlist2
                            u = i < j ? i : j
                            v = i < j ? j : i
                            push!(rows[u], v)
                        end
                        j = nextp[j]
                    end
                end
                i = nextp[i]
            end
        end
    elseif method == :bruteforce
        rlist2 = (r_verlet + master_nl.skin)^2
        @inbounds for i in 1:N-1
            for j in i+1:N
                r2 = _squared_distance_min_image(R, i, j, box)
                if r2 <= rlist2
                    push!(rows[i], j)
                end
            end
        end
    elseif method == :all_pairs
        @inbounds for i in 1:N-1
            for j in i+1:N
                push!(rows[i], j)
            end
        end
    else
        error("Unknown neighborlist method: $method")
    end

    # Build CSR from row buffers
    rowptr = master_nl.rowptr
    if length(rowptr) != N + 1
        resize!(rowptr, N + 1)
    end
    rowptr[1] = 1
    @inbounds for i in 1:N
        rowptr[i+1] = rowptr[i] + length(rows[i])
    end
    M = rowptr[end] - 1

    if length(master_nl.colind) < M
        newlen = max(M, max(1, length(master_nl.colind)) * 2)
        resize!(master_nl.colind, newlen)
    end

    colind = master_nl.colind
    idx = 1
    @inbounds for i in 1:N
        li = rows[i]
        leni = length(li)
        if leni > 0
            copyto!(colind, idx, li, 1, leni)
            idx += leni
        end
    end

    return master_nl
end
# function build_master_neighborlist!(master_nl::MasterNeighborCSRList, R::AbstractVector, box::CubicBox;
#     r_verlet::Real=0.0, method::Symbol=:cells, grid::Union{Nothing,CellGrid}=nothing)

#     @assert length(R[1]) == 3 "d=3 only"
#     N = length(R)

#     # Buffer for accepted pairs (u<i<j, v=j); reuse across steps to avoid reallocs
#     pairs = master_nl.pairs
#     empty!(pairs)

#     # Single-sweep collection into per-row buffers to avoid double distance work (reuse scratch)
#     rows = master_nl.rows
#     if length(rows) != N
#         resize!(rows, N)
#         for i in 1:N
#             empty!(rows[i])
#         end
#     else
#         for i in 1:N
#             empty!(rows[i])
#         end
#     end

#     if method == :cells
#         rlist2 = (r_verlet + master_nl.skin)^2
#         rlist = sqrt(rlist2)
#         g = if grid === nothing || grid.cell_size < rlist || grid.L != box.L
#             build_cellgrid(R, box; cell_size=rlist)
#         else
#             rebin!(grid, R, box)
#         end

#         nx, ny, nz = Int.(g.dims)
#         heads = g.heads; nextp = g.next
#         @inbounds for cz in 1:nz, cy in 1:ny, cx in 1:nx
#             c_idx = ((cz - 1) * ny + (cy - 1)) * nx + cx
#             neighs = StaticArrays.MVector{13,Int}(undef); ncount = 0
#             for dz in -1:1, dy in -1:1, dx in 0:1
#                 if (dx == 0 && dy < 0) || (dx == 0 && dy == 0 && dz <= 0)
#                     continue
#                 end
#                 ncx = mod1(cx + dx, nx); ncy = mod1(cy + dy, ny); ncz = mod1(cz + dz, nz)
#                 cc_idx = ((ncz - 1) * ny + (ncy - 1)) * nx + ncx
#                 if cc_idx == c_idx; continue; end
#                 seen = false
#                 @inbounds for t in 1:ncount
#                     if neighs[t] == cc_idx; seen = true; break; end
#                 end
#                 if !seen
#                     ncount += 1; neighs[ncount] = cc_idx
#                 end
#             end

#             i = heads[c_idx]
#             while i != 0
#                 j = nextp[i]
#                 while j != 0
#                     r2 = _squared_distance_min_image(R, i, j, box)
#                     if r2 <= rlist2
#                         u = i < j ? i : j
#                         v = i < j ? j : i
#                         push!(pairs, (u, v))
#                     end
#                     j = nextp[j]
#                 end
#                 @inbounds for t in 1:ncount
#                     cc_idx = neighs[t]
#                     j = heads[cc_idx]
#                     while j != 0
#                         r2 = _squared_distance_min_image(R, i, j, box)
#                         if r2 <= rlist2
#                             u = i < j ? i : j
#                             v = i < j ? j : i
#                             push!(pairs, (u, v))
#                         end
#                         j = nextp[j]
#                     end
#                 end
#                 i = nextp[i]
#             end
#         end
#     elseif method == :bruteforce
#         rlist2 = (r_verlet + master_nl.skin)^2
#         @inbounds for i in 1:N-1
#             for j in i+1:N
#                 r2 = _squared_distance_min_image(R, i, j, box)
#                 if r2 <= rlist2
#                     push!(rows[i], j)
#                 end
#             end
#         end
#     elseif method == :all_pairs
#         @inbounds for i in 1:N-1
#             for j in i+1:N
#                 push!(rows[i], j)
#             end
#         end
#     else
#         error("Unknown neighborlist method: $method")
#     end

#     # Build CSR from pairs: counts → rowptr → fill colind
#     counts = zeros(Int, N)
#     @inbounds for (u, _) in pairs
#         counts[u] += 1
#     end
#     rowptr = Vector{Int}(undef, N + 1)
#     rowptr[1] = 1
#     @inbounds for i in 1:N
#         rowptr[i+1] = rowptr[i] + counts[i]
#     end
#     M = rowptr[end] - 1
#     colind = Vector{Int}(undef, M)
#     writeptr = copy(rowptr)
#     @inbounds for (u, v) in pairs
#         k = writeptr[u]
#         colind[k] = v
#         writeptr[u] = k + 1
#     end

#     master_nl.rowptr = rowptr
#     master_nl.colind = colind
#     return master_nl
# end



"""
    NeighborList{IT<:Integer,T<:Real}

Compressed-sparse-row (CSR) neighbor list with snapshot of reference positions.

Fields
- `cutoff::T`: physical cutoff distance (`rcut`).
- `skin::T`: extra buffer distance beyond cutoff.
- `pairs::Vector{IT}`: concatenated neighbor indices (CSR-style).
- `offsets::Vector{IT}`: length `N+1`; neighbors of particle `i` are stored in `pairs[offsets[i]:(offsets[i+1]-1)]`.
- `ref_positions::Matrix{T}`: snapshot of positions at last rebuild, used for displacement checks.

Use `build_neighborlist` to construct, and `maybe_rebuild!` to keep it updated.
"""
mutable struct NeighborList{IT,T, Dims}
    cutoff::T
    skin::T
    pairs::Vector{IT}         # concatenated neighbor indices (1-based)
    offsets::Vector{IT}       # length N+1; neighbors of i are pairs[offsets[i]:(offsets[i+1]-1)]
    ref_positions::Vector{SVector{Dims,T}}  # N copy at last (re)build
end

"""
    build_neighborlist(R, box; cutoff, skin=0.3) -> NeighborList

Construct a neighbor list for positions `R` inside a `CubicBox`.
Pairs are included if they are within `cutoff + skin` under minimum-image distance.

This is an **O(N²)** operation but performed infrequently. The returned `NeighborList`
uses a CSR (compressed sparse row) layout for compact storage and fast iteration.
"""
function build_neighborlist(R::Vector{SVector{Dims,T}}, box::CubicBox{T}; cutoff::T, skin::T=T(0.3)) where {Dims,T}
    N = length(R)
    IT = Int
    rlist2 = (cutoff + skin)^2
    neighs = [Vector{IT}() for _ in 1:N]
    Δ = zeros(T, Dims)
    @inbounds for i in 1:N-1
        ri = R[i]
        for j in i+1:N
            rj = R[j]
            Δ .= ri - rj
            minimum_image!(Δ, box)
            r2 = dot(Δ, Δ)
            if r2 <= rlist2
                push!(neighs[i], IT(j))
                push!(neighs[j], IT(i))
            end
        end
    end
    offsets = Vector{IT}(undef, N+1)
    offsets[1] = IT(1)
    @inbounds for i in 1:N
        offsets[i+1] = offsets[i] + IT(length(neighs[i]))
    end
    pairs = Vector{IT}(undef, Int(offsets[end])-1)
    @inbounds for i in 1:N
        start = offsets[i]
        stop  = offsets[i+1]-1
        if stop >= start
            pairs[start:stop] = neighs[i]
        end
    end
    NeighborList{IT,T, Dims}(cutoff, skin, pairs, offsets, copy(R))
end

"""
    max_displacement_since_build(nlist, R, box) -> Float64

Return the maximum particle displacement (with minimum-image convention) since `nlist` was built.
Used internally by `maybe_rebuild!`.
"""
function max_displacement_since_build(nlist::NeighborList{IT,T, Dims}, R::Vector{SVector{Dims,T}}, box::CubicBox{T}) where {IT,T,Dims}
    N = length(R)
    @assert length(nlist.ref_positions) == N
    Δ = zeros(T, Dims)
    maxdisp = 0.0
    @inbounds for i in 1:N
        ri = R[i]
        r0 = nlist.ref_positions[i]
        Δ .= ri - r0
        minimum_image!(Δ, box)
        d = maximum(abs, Δ)
        if d > maxdisp
            maxdisp = d
        end
    end
    return maxdisp
end

"""
    maybe_rebuild!(nlist, R, box) -> Bool

Check whether the neighbor list `nlist` is still valid given new positions `R`.
If the maximum displacement since the last build is greater than `skin/2`, rebuild the
list and return `true`. Otherwise leave it unchanged and return `false`.
"""
function maybe_rebuild!(nlist::NeighborList{IT,T, Dims}, R::Vector{SVector{Dims,T}}, box::CubicBox{T}) where {IT,T,Dims}
    maxdisp = max_displacement_since_build(nlist, R, box)
    if maxdisp > nlist.skin / 2
        newnl = build_neighborlist(R, box; cutoff=nlist.cutoff, skin=nlist.skin)
        nlist.pairs = newnl.pairs
        nlist.offsets = newnl.offsets
        nlist.ref_positions = newnl.ref_positions
        return true
    else
        return false
    end
end

"""
    lj_forces(R, box, nlist; ϵ=1.0, σ=1.0, rcut=NaN, shift=false, return_potential=false)

Compute Lennard–Jones forces using a prebuilt [`NeighborList`](@ref).
This is more efficient than the brute-force O(N²) kernel, scaling ~O(N) for fixed density.

Arguments
- `R::AbstractMatrix`: positions (N×D).
- `box::CubicBox`: periodic cubic box.
- `nlist::NeighborList`: neighbor list built with at least the requested cutoff.

Keywords
- `ϵ`: LJ well depth (default 1.0).
- `σ`: LJ length scale (default 1.0).
- `rcut`: optional cutoff to use for the force calculation (defaults to `nlist.cutoff`).
- `shift::Bool`: if true, subtract potential at cutoff to make it continuous.
- `return_potential::Bool`: if true, also return total potential energy.

Returns
- `F::Matrix`: force on each particle (same shape as `R`).
- `U::Float64` (if `return_potential=true`).

Example
```@example
box = CubicBox(8.0)
R = randn(10,3)
wrap_positions!(R, box)
nlist = build_neighborlist(R, box; cutoff=2.5, skin=0.4)
F, U = lj_forces(R, box, nlist; rcut=2.5, return_potential=true)
```
"""
function lj_forces(R::Vector{SVector{Dims,T}}, box::CubicBox{T}, nlist::NeighborList{IT,T, Dims};
                   ϵ::Real=1.0, σ::Real=1.0, rcut::Real=NaN,
                   shift::Bool=false, return_potential::Bool=false) where {Dims,T,IT}
    N = length(R)
    F = [zero(SVector{Dims,T}) for _ in 1:N]
    U = 0.0
    σ2 = float(σ)^2
    cutoff = isnan(rcut) ? nlist.cutoff : float(rcut)
    rcut2 = cutoff^2

    # Energy shift at cutoff (no smoothing of forces)
    Uc = 0.0
    if shift && isfinite(cutoff)
        s2c = σ2 / rcut2
        s6c = s2c^3
        Uc  = 4*float(ϵ)*(s6c^2 - s6c)
    end

    Δ = zeros(T, Dims)
    @inbounds for i in 1:N
        ri = R[i]
        for idx in nlist.offsets[i]:(nlist.offsets[i+1]-1)
            j = nlist.pairs[idx]
            # Apply each symmetric pair once
            if j > i
                rj = R[j]
                Δ .= ri - rj
                minimum_image!(Δ, box)
                r2 = dot(Δ, Δ)
                if r2 == 0.0
                    continue
                end
                if r2 <= rcut2
                    invr2 = 1 / r2
                    s2 = σ2 * invr2
                    s6 = s2^3
                    fr_over_r = 24*float(ϵ)*(2*s6^2 - s6) * invr2
                    fvec = fr_over_r .* Δ
                    F[i] += fvec
                    F[j] -= fvec
                    if return_potential
                        U += 4*float(ϵ)*(s6^2 - s6) - Uc
                    end
                end
            end
        end
    end
    return return_potential ? (F, U) : F
end




function lj_forces(positions::AbstractMatrix, box::CubicBox;
                   ϵ::Real=1.0, σ::Real=1.0, rcut::Real=Inf,
                   shift::Bool=false, return_potential::Bool=false)
    N, D = size(positions)
    F = zeros(Float64, N, D)
    U = 0.0

    σ2 = float(σ)^2
    rcut2 = float(rcut)^2

    Uc = 0.0
    if shift && isfinite(rcut)
        s2c = σ2 / rcut2
        s6c = s2c^3
        Uc  = 4*float(ϵ)*(s6c^2 - s6c)
    end

    Δ = zeros(Float64, D)

    @inbounds for i in 1:N-1
        ri = @view positions[i, :]
        for j in i+1:N
            rj = @view positions[j, :]
            Δ .= ri .- rj
            minimum_image!(Δ, box)
            r2 = dot(Δ, Δ)
            if r2 == 0.0
                continue
            end
            if r2 <= rcut2
                invr2 = 1 / r2
                s2 = σ2 * invr2
                s6 = s2^3
                fr_over_r = 24*float(ϵ)*(2*s6^2 - s6) * invr2
                for k in 1:D
                    f = fr_over_r * Δ[k]
                    F[i,k] += f
                    F[j,k] -= f
                end
                if return_potential
                    U += 4*float(ϵ)*(s6^2 - s6) - Uc
                end
            end
        end
    end
    return return_potential ? (F, U) : F
end

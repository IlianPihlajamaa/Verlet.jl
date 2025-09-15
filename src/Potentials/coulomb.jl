using LinearAlgebra, StaticArrays

mutable struct Coulomb{IntT<:Integer, T<:Core.AbstractPotentialPair, T_Float} <: Core.AbstractPairPotential
    params::PairTable{T}
    exclusions::Vector{Tuple{IntT, IntT}}
    neighbors::Neighbors.PotentialNeighborList{T}
    skin::T_Float
end

function Coulomb(params::PairTable{T}, exclusions, skin) where T
    neighbors = StructArray{Neighbors.NeighborPair{T, Core.T_Int}}(undef, 0)
    return Coulomb(params, exclusions, neighbors, skin)
end

function Core.compute_forces!(pot::Coulomb, sys::Core.System)
    for pair_info in pot.neighbors
        i = pair_info.i
        j = pair_info.j
        p = pair_info.pair

        Δ = sys.positions[i] - sys.positions[j]
        Δ = Core.minimum_image(Δ, sys.box)
        r2 = dot(Δ, Δ)

        if r2 < p.rc^2
            q1q2 = p.q1q2
            invr2 = 1 / r2
            r = sqrt(r2)
            invr = 1 / r
            fr_over_r = q1q2 * invr2 * invr # k=1
            fvec = fr_over_r .* Δ
            sys.forces[i] += fvec
            sys.forces[j] -= fvec
        end
    end
end

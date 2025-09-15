using LinearAlgebra, StaticArrays



mutable struct LennardJones{IntT<:Integer, T<:Core.AbstractPotentialPair, T_Float} <: Core.AbstractPairPotential
    params::PairTable{T}
    exclusions::Vector{Tuple{IntT, IntT}}
    neighbors::Neighbors.PotentialNeighborList{T}
    skin::T_Float
end

function LennardJones(params::PairTable{T}, exclusions, skin) where T
    neighbors = StructArray{Neighbors.NeighborPair{T, Core.T_Int}}(undef, 0)
    return LennardJones(params, exclusions, neighbors, skin)
end

function Core.compute_forces!(pot::LennardJones, sys::Core.System)
    for pair_info in pot.neighbors
        i = pair_info.i
        j = pair_info.j
        p = pair_info.pair

        Δ = sys.positions[i] - sys.positions[j]
        Δ = Core.minimum_image(Δ, sys.box)
        r2 = dot(Δ, Δ)

        if r2 < p.rc^2
            ϵ = p.ε
            σ = p.σ
            σ2 = σ^2
            invr2 = 1 / r2
            s2 = σ2 * invr2
            s6 = s2^3
            fr_over_r = 24 * ϵ * (2 * s6^2 - s6) * invr2
            fvec = fr_over_r .* Δ
            sys.forces[i] += fvec
            sys.forces[j] -= fvec
        end
    end
end

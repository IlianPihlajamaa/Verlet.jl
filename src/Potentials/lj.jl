import ..Core: compute_forces!, compute_potential_energy

mutable struct LennardJones{IntT<:Integer, T<:AbstractPotentialPair, T_Float} <: AbstractPairPotential
    params::PairTable{T}
    exclusions::Vector{Tuple{IntT, IntT}}
    neighborlist::PotentialNeighborList{T, IntT}
    skin::T_Float
end

function LennardJones(params::PairTable{T}, exclusions, skin) where T
    neighborlist = PotentialNeighborList(eltype(params.table))
    return LennardJones(params, exclusions, neighborlist, skin)
end



function compute_forces!(pot::LennardJones, sys::System)
    for pair_info in pot.neighborlist.neighbors
        i = pair_info.i
        j = pair_info.j
        p = pair_info.pair

        Δ = sys.positions[i] - sys.positions[j]
        Δ = minimum_image(Δ, sys.box)
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

function compute_potential_energy(pot::LennardJones, sys::System)
    energy = zero(eltype(sys.positions[1]))
    for pair_info in pot.neighborlist.neighbors
        i = pair_info.i
        j = pair_info.j
        params = pair_info.pair

        Δ = sys.positions[i] - sys.positions[j]
        Δ = minimum_image(Δ, sys.box)
        r2 = dot(Δ, Δ)

        if r2 < params.rc^2
            invr2 = 1 / r2
            s2 = (params.σ^2) * invr2
            s6 = s2^3
            energy += 4 * params.ε * (s6^2 - s6)
        end
    end
    return energy
end

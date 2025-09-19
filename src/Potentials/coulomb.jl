import ..Core: compute_forces!, compute_potential_energy
mutable struct Coulomb{IntT<:Integer, T<:AbstractPotentialPair, T_Float} <: AbstractPairPotential
    params::PairTable{T}
    exclusions::Vector{Tuple{IntT, IntT}}
    neighborlist::PotentialNeighborList{T}
    skin::T_Float
end

function Coulomb(params::PairTable{T}, exclusions, skin) where T
    neighborlist = PotentialNeighborList(eltype(params.table))
    return Coulomb(params, exclusions, neighborlist, skin)
end

function compute_forces!(pot::Coulomb, sys::System)
    for pair_info in pot.neighborlist.neighbors
        i = pair_info.i
        j = pair_info.j
        p = pair_info.pair

        Δ = sys.positions[i] - sys.positions[j]
        Δ = minimum_image(Δ, sys.box)
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

function compute_potential_energy(pot::Coulomb, sys::System)
    energy = zero(eltype(sys.positions[1]))
    for pair_info in pot.neighborlist.neighbors
        i = pair_info.i
        j = pair_info.j
        params = pair_info.pair

        Δ = sys.positions[i] - sys.positions[j]
        Δ = minimum_image(Δ, sys.box)
        r2 = dot(Δ, Δ)

        if r2 < params.rc^2 && r2 > 0
            energy += params.q1q2 / sqrt(r2)
        end
    end
    return energy
end

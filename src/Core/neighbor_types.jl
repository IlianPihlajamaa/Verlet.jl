using StructArrays

struct NeighborPair{F<:AbstractPotentialPair, IntT<:Integer}
    i::IntT
    j::IntT
    pair::F
end

const PotentialNeighborList{F} = StructArray{NeighborPair{F, T_int}} where {F<:AbstractPotentialPair}

struct MasterNeighborEntry
    i::T_int
    j::T_int
    r2::T_Float
end

mutable struct MasterNeighborList{T<:Number}
    skin::T
    entries::Vector{MasterNeighborEntry}
    # Other fields like reference positions could be added here if needed
end

function MasterNeighborList(skin::T; sizehint=1000) where T
    entries = MasterNeighborEntry[]
    sizehint!(entries, sizehint)
    return MasterNeighborList(skin, entries)
end

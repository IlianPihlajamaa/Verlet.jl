struct NeighborPair{F<:AbstractPotentialPair, IntT<:Integer}
    i::IntT
    j::IntT
    pair::F
end

struct PotentialNeighborList{F, T_Int} 
    neighbors::Vector{NeighborPair{F, T_Int}}
    PotentialNeighborList(::Type{F}) where F = new{F, T_Int}(NeighborPair{F, T_Int}[])
end

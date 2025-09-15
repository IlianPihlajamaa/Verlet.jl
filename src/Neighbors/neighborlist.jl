struct NeighborPair{F<:AbstractPotentialPair, IntT<:Integer}
    i::IntT
    j::IntT
    pair::F
end

struct PotentialNeighborList{F, T_Int} 
    neighbors::Vector{NeighborPair{F, T_Int}}
    PotentialNeighborList(::Type{F}) where F = new{F, T_Int}(NeighborPair{F, T_Int}[])
end

struct MasterNeighborEntry{T_Float, T_Int<:Integer}
    i::T_Int
    j::T_Int
    r2::T_Float
end

mutable struct MasterNeighborList{T<:Number, T_Int<:Integer} <: AbstractNeighborList
    skin::T
    entries::Vector{MasterNeighborEntry{T, T_Int}}
    # Other fields like reference positions could be added here if needed
end

function MasterNeighborList(skin::T; sizehint=1000) where T
    entries = MasterNeighborEntry{T, Int}[]  # Specify the type parameters
    sizehint!(entries, sizehint)
    return MasterNeighborList(skin, entries)
end

# CSR-style master neighbor list (topology only, no r2)
mutable struct MasterNeighborCSRList{T<:Number, IT<:Integer} <: AbstractNeighborList
    skin::T
    rowptr::Vector{IT}                 # length N+1; neighbors of i in rowptr[i]:(rowptr[i+1]-1)
    colind::Vector{IT}                 # concatenation of j (> i) indices
    counts::Vector{IT}                 # scratch counts per row (preallocated, zeroed each build)
    writeptr::Vector{IT}               # scratch write positions (same length as rowptr)
    rows::Vector{Vector{IT}}           # legacy scratch; kept for compatibility
    pairs::Vector{NTuple{2,Int}}       # legacy scratch; kept for compatibility
end

function MasterNeighborCSRList(skin::T; N::Integer=0, M_hint::Integer=0) where {T}
    N = max(0, Int(N))
    rowptr = fill(Int(1), N + 1)
    colind = Vector{Int}(undef, max(0, Int(M_hint)))
    counts = Vector{Int}(undef, N)
    fill!(counts, 0)
    writeptr = similar(rowptr)
    rows = [Int[] for _ in 1:N]
    pairs = NTuple{2,Int}[]
    return MasterNeighborCSRList{T, Int}(skin, rowptr, colind, counts, writeptr, rows, pairs)
end

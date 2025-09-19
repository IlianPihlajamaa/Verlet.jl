module Neighbors

using ..Core
using LinearAlgebra, StaticArrays, StructArrays
import ..Core: rebuild_neighbors!, prepare_neighbors!, System

include("neighborlist.jl")
include("cellgrid.jl")
include("skin_master_nl.jl")
include("forcefields.jl")

export build_master_neighborlist!, build_cellgrid, rebin!, ForceField, build_all_neighbors!, compute_all_forces!
export PotentialNeighborList, MasterNeighborList, brute_force_pairs

function rebuild_neighbors!(system::System, master::MasterNeighborList; kwargs...)
    ff = system.forcefield
    ff === nothing && throw(ArgumentError("System has no forcefield attached; cannot rebuild neighbor lists"))
    build_all_neighbors!(master, ff, system; kwargs...)
    return master
end

end

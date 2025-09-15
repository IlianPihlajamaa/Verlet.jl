using Test
using Verlet

using LinearAlgebra, StaticArrays

@testset "Core" begin
    include("test_system.jl")
    include("test_core.jl")
end

@testset "Neighbors" begin
    # include("test_neighborlist.jl") # Old API
    include("test_cellgrid.jl")
    include("test_forcefields.jl")
end

@testset "Thermostats" begin
    include("test_thermostats.jl")
end

@testset "Constraints" begin
    include("test_constraints.jl")
	include("test_cbaoab.jl")
end

@testset "Potentials" begin
    include("test_bonded.jl")
end

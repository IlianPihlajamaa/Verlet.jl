module VerletTests

using Test
using Verlet

using LinearAlgebra, StaticArrays

# Core unit tests
@testset "Core" begin
    include("test_system.jl")
    include("test_core.jl")
end

# Neighborlist and cellgrid acceptance tests
@testset "Neighbors" begin
    # include("test_neighborlist.jl") # Old API
    include("test_cellgrid.jl")
    include("test_forcefields.jl")
end

# Thermostat / NVT tests
@testset "Thermostats" begin
    include("test_thermostats.jl")
end

# Constraints tests
@testset "Constraints" begin
    include("test_constraints.jl")
	include("test_cbaoab.jl")
end

# Bonded potential tests
@testset "Potentials" begin
    include("test_bonded.jl")
end

end

using Verlet
using StaticArrays
using Test

@testset "System struct" begin
    # 3D system
    positions = [@SVector randn(3) for _ in 1:10]
    velocities = [@SVector randn(3) for _ in 1:10]
    forces = [@SVector zeros(3) for _ in 1:10]
    masses = ones(10)
    box = CubicBox(10.0)
    types = ones(Int, 10)
    type_names = Dict(1 => :A)

    # Valid constructor
    sys = System(positions, velocities, forces, masses, box, types, type_names)
    @test natoms(sys) == 10
    @test natomtypes(sys) == 1

    # Inconsistent sizes
    @test_throws AssertionError System(positions, velocities[1:9], forces, masses, box, types, type_names)
    @test_throws AssertionError System(positions, velocities, forces[1:9], masses, box, types, type_names)
    @test_throws AssertionError System(positions, velocities, forces, masses[1:9], box, types, type_names)
    @test_throws AssertionError System(positions, velocities, forces, masses, box, types[1:9], type_names)

    # 2D system
    positions2d = [@SVector randn(2) for _ in 1:5]
    velocities2d = [@SVector randn(2) for _ in 1:5]
    forces2d = [@SVector zeros(2) for _ in 1:5]
    masses2d = ones(5)
    box2d = CubicBox(5.0)
    types2d = ones(Int, 5)
    type_names2d = Dict(1 => :A)

    sys2d = System(positions2d, velocities2d, forces2d, masses2d, box2d, types2d, type_names2d)
    @test natoms(sys2d) == 5
    @test natomtypes(sys2d) == 1
end

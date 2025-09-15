using Test
using Verlet
using StaticArrays
using LinearAlgebra, StaticArrays

const Dims = 3

@testset "Core functionality" begin
    @testset "System construction" begin
        positions = [SVector{Dims,Float64}(0.0,0.0,0.0)]
        velocities = [SVector{Dims,Float64}(1.0,0.0,0.0)]
        forces = [SVector{Dims,Float64}(0.0,0.0,0.0)]
        masses = [1.0]
        box = CubicBox(10.0)
        types = [1]
        type_names = Dict(1 => :A)
        sys = System(positions, velocities, forces, masses, box, types, type_names)
        @test length(sys.positions) == 1
        @test length(sys.velocities) == 1
        @test sys.masses == [1.0]
    end
    @testset "Velocity Verlet" begin
        @testset "Velocity Verlet free particle" begin
            forces_func(r) = fill(SVector{Dims,Float64}(0.0,0.0,0.0), length(r))
            positions = [SVector{Dims,Float64}(0.0,0.0,0.0)]
            velocities = [SVector{Dims,Float64}(1.0,0.0,0.0)]
            forces = [SVector{Dims,Float64}(0.0,0.0,0.0)]
            masses = [1.0]
            box = CubicBox(10.0)
            types = [1]
            type_names = Dict(1 => :A)
            sys = System(positions, velocities, forces, masses, box, types, type_names)
            velocity_verlet!(sys, forces_func, 0.1)
            @test isapprox(sys.positions[1][1], 0.1; atol=1e-12)
            @test isapprox(sys.velocities[1][1], 1.0; atol=1e-12)
        end

        @testset "Velocity Verlet harmonic oscillator" begin
            forces_func(r) = map(x -> -x, r)
            positions = [SVector{Dims,Float64}(1.0,0.0,0.0)]
            velocities = [SVector{Dims,Float64}(0.0,0.0,0.0)]
            forces = [SVector{Dims,Float64}(0.0,0.0,0.0)]
            masses = [1.0]
            box = CubicBox(10.0)
            types = [1]
            type_names = Dict(1 => :A)
            sys = System(positions, velocities, forces, masses, box, types, type_names)
            velocity_verlet!(sys, forces_func, 0.1)
            @test sys.positions[1][1] < 1.0   # should move left
        end
    end

    @testset "kinetic energy" begin
        # v = (3,4,0) => v^2 = 25; KE = 0.5 * m * v^2
        positions = [SVector{Dims,Float64}(0.0,0.0,0.0)]
        velocities = [SVector{Dims,Float64}(3.0,4.0,0.0)]
        forces = [SVector{Dims,Float64}(0.0,0.0,0.0)]
        masses = [2.0]
        box = CubicBox(10.0)
        types = [1]
        type_names = Dict(1 => :A)
        sys = System(positions, velocities, forces, masses, box, types, type_names)
        @test isapprox(kinetic_energy(sys), 0.5 * 2.0 * 25.0; atol=1e-12)
    end

    @testset "potential energy conventions" begin
        # A force that supports keyword return_potential
        function f_kw(r; return_potential::Bool=false)
            F = map(x -> -x, r)
            if return_potential

                U = 0.5 * sum(([sum(abs2.(ri)) for ri in r]))  # harmonic potential: 0.5 * Σ r^2
                return (F, U)
            else
                return F
            end
        end
        positions = [SVector{Dims,Float64}(1.0,0.0,0.0)]
        velocities = [SVector{Dims,Float64}(0.0,0.0,0.0)]
        forces = [SVector{Dims,Float64}(0.0,0.0,0.0)]
        masses = [1.0]
        box = CubicBox(10.0)
        types = [1]
        type_names = Dict(1 => :A)
        sys = System(positions, velocities, forces, masses, box, types, type_names)
        @test isapprox(Verlet.Core.potential_energy(sys, f_kw), 0.5; atol=1e-12)

        # A force that does NOT provide potential should error
        f_plain(r) = map(x -> -x, r)
        @test_throws ErrorException Verlet.Core.potential_energy(sys, f_plain)
    end

    @testset "CubicBox + minimum_image!" begin
        box = CubicBox(10.0)
        Δ = [ 6.0, -6.0,  0.1 ]  # in 3D
        Δ = minimum_image(Δ, box)
        @test Δ[1] ==  6.0 - 10.0    # -> -4.0
        @test Δ[2] == -6.0 + 10.0    # ->  4.0
        @test isapprox(Δ[3], 0.1; atol=1e-12)
    end
end
using Test
using Verlet
using StaticArrays
using LinearAlgebra, StaticArrays

const Dims = 3

@testset "Core functionality" begin
    @testset "ParticleSystem construction" begin
        ps = ParticleSystem{Dims,Float64}([SVector{Dims,Float64}(0.0,0.0,0.0)], [SVector{Dims,Float64}(1.0,0.0,0.0)], [1.0])
        @test length(ps.positions) == 1
        @test length(ps.velocities) == 1
        @test ps.masses == [1.0]
    end
    @testset "Velocity Verlet" begin
        @testset "Velocity Verlet free particle" begin
            forces(r) = fill(SVector{Dims,Float64}(0.0,0.0,0.0), length(r))
            ps = ParticleSystem{Dims,Float64}([SVector{Dims,Float64}(0.0,0.0,0.0)], [SVector{Dims,Float64}(1.0,0.0,0.0)], [1.0])
            velocity_verlet!(ps, forces, 0.1)
            @test isapprox(ps.positions[1][1], 0.1; atol=1e-12)
            @test isapprox(ps.velocities[1][1], 1.0; atol=1e-12)
        end

        @testset "Velocity Verlet harmonic oscillator" begin
            forces(r) = map(x -> -x, r)
            ps = ParticleSystem{Dims,Float64}([SVector{Dims,Float64}(1.0,0.0,0.0)], [SVector{Dims,Float64}(0.0,0.0,0.0)], [1.0])
            velocity_verlet!(ps, forces, 0.1)
            @test ps.positions[1][1] < 1.0   # should move left
        end
    end

    @testset "kinetic energy" begin
        # v = (3,4,0) => v^2 = 25; KE = 0.5 * m * v^2
        ps = ParticleSystem{Dims,Float64}([SVector{Dims,Float64}(0.0,0.0,0.0)], [SVector{Dims,Float64}(3.0,4.0,0.0)], [2.0])
        @test isapprox(kinetic_energy(ps), 0.5 * 2.0 * 25.0; atol=1e-12)
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
        ps = ParticleSystem{Dims,Float64}([SVector{Dims,Float64}(1.0,0.0,0.0)], [SVector{Dims,Float64}(0.0,0.0,0.0)], [1.0])
        @test isapprox(potential_energy(ps, f_kw), 0.5; atol=1e-12)

        # A force that does NOT provide potential should error
        f_plain(r) = map(x -> -x, r)
        @test_throws ErrorException potential_energy(ps, f_plain)
    end

    @testset "CubicBox + minimum_image!" begin
        box = CubicBox(10.0)
        Δ = [ 6.0, -6.0,  0.1 ]  # in 3D
        Δ = minimum_image(Δ, box)
        @test Δ[1] ==  6.0 - 10.0    # -> -4.0
        @test Δ[2] == -6.0 + 10.0    # ->  4.0
        @test isapprox(Δ[3], 0.1; atol=1e-12)
    end

    # The following tests for LJ and neighbor lists will need to be refactored for SVector-based positions.
    # ...existing code...
end
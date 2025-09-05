using Test
using Test
using Verlet

@testset "ParticleSystem construction" begin
    ps = ParticleSystem([0.0 0.0], [1.0 0.0], [1.0])
    @test size(ps.positions) == (1,2)
    @test size(ps.velocities) == (1,2)
    @test ps.masses == [1.0]
end

@testset "Velocity Verlet free particle" begin
    forces(r) = zeros(size(r))
    ps = ParticleSystem([0.0 0.0], [1.0 0.0], [1.0])
    velocity_verlet!(ps, forces, 0.1)
    @test isapprox(ps.positions[1,1], 0.1; atol=1e-12)
    @test isapprox(ps.velocities[1,1], 1.0; atol=1e-12)
end

@testset "Velocity Verlet harmonic oscillator" begin
    forces(r) = -r
    ps = ParticleSystem([1.0 0.0], [0.0 0.0], [1.0])
    velocity_verlet!(ps, forces, 0.1)
    @test ps.positions[1,1] < 1.0   # should move left
end

@testset "kinetic energy" begin
    # v = (3,4) => v^2 = 25; KE = 0.5 * m * v^2
    ps = ParticleSystem([0.0 0.0], [3.0 4.0], [2.0])
    @test isapprox(kinetic_energy(ps), 0.5 * 2.0 * 25.0; atol=1e-12)
end

@testset "potential energy conventions" begin
    # A force that supports keyword return_potential
    function f_kw(r; return_potential::Bool=false)
        F = -r
        if return_potential
            U = 0.5 * sum(r.^2)  # harmonic potential: 0.5 * Σ r^2
            return (F, U)
        else
            return F
        end
    end
    ps = ParticleSystem([1.0 0.0], [0.0 0.0], [1.0])
    @test isapprox(potential_energy(ps, f_kw), 0.5; atol=1e-12)

    # A force that does NOT provide potential should error
    f_plain(r) = -r
    @test_throws ErrorException potential_energy(ps, f_plain)
end
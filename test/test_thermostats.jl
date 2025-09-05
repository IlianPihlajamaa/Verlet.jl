using Test
using Random
using Verlet

@testset "thermostats.jl" begin
    @testset "instantaneous_temperature & degrees_of_freedom" begin
        N, D = 10, 3
        ps = ParticleSystem(zeros(N,D), zeros(N,D), ones(N))
        # set velocities so KE = 0.5 * dof * kB*T with kB=1, T=2
        Tt, kB = 2.0, 1.0
        dof = N*D
        KE_target = 0.5 * dof * kB * Tt
        # put all energy in v[:,1]
        vmag = sqrt(2*KE_target / sum(ps.masses))
        ps.velocities[:,1] .= vmag
        @test degrees_of_freedom(ps) == dof
        @test isapprox(instantaneous_temperature(ps; kB=kB), Tt; rtol=1e-12, atol=0)
    end

    @testset "velocity_rescale! matches target T" begin
        N, D = 8, 3
        ps = ParticleSystem(zeros(N,D), randn(N,D), ones(N))
        Tt = 3.0
        velocity_rescale!(ps, Tt; kB=1.0)
        @test isapprox(instantaneous_temperature(ps; kB=1.0), Tt; rtol=1e-12, atol=0)
    end

    @testset "BAOAB: γ=0 reduces to VV trajectory (single step)" begin
        # Compare one BAOAB step with γ=0 against velocity_verlet! for a harmonic oscillator.
        using LinearAlgebra
        N, D = 4, 3
        k = 1.0
        ho_forces(R) = -k .* R
        ps1 = ParticleSystem(randn(N,D), randn(N,D), ones(N))
        ps2 = ParticleSystem(copy(ps1.positions), copy(ps1.velocities), copy(ps1.masses))
        dt = 0.01
        # BAOAB with γ=0 and a fixed RNG (unused)
        rng = MersenneTwister(42)
        langevin_baoab!(ps1, ho_forces, dt; γ=0.0, T=0.0, kB=1.0, rng=rng)
        velocity_verlet!(ps2, ho_forces, dt)
        @test isapprox(ps1.positions, ps2.positions; rtol=1e-12, atol=1e-12)
        @test isapprox(ps1.velocities, ps2.velocities; rtol=1e-12, atol=1e-12)
    end

    @testset "BAOAB: temperature control converges to target (statistical)" begin
        # Stochastic test with tolerance: run several steps and check mean T near target.
        N, D = 64, 3
        ps = ParticleSystem(zeros(N,D), randn(N,D), ones(N))
        Tt, γ, dt = 1.5, 1.0, 0.005
        rng = MersenneTwister(123)
        steps = 2000
        for _ in 1:steps
            # harmonic bath keeps system bound; any stable force will do
            ho_forces(R) = -0.1 .* R
            langevin_baoab!(ps, ho_forces, dt; γ=γ, T=Tt, kB=1.0, rng=rng)
        end
        Tinst = instantaneous_temperature(ps; kB=1.0)
        # Allow ~5% tolerance; stochastic fluctuations decrease with N and steps
        @test isapprox(Tinst, Tt; rtol=0.05)
    end

    @testset "BAOAB reproducibility with fixed RNG" begin
        N, D = 16, 3
        psA = ParticleSystem(randn(N,D), randn(N,D), ones(N))
        psB = ParticleSystem(copy(psA.positions), copy(psA.velocities), copy(psA.masses))
        Tt, γ, dt = 2.0, 0.7, 0.01
        rngA = MersenneTwister(2025)
        rngB = MersenneTwister(2025)
        steps = 50
        for _ in 1:steps
            f(R) = -0.2 .* R
            langevin_baoab!(psA, f, dt; γ=γ, T=Tt, kB=1.0, rng=rngA)
            langevin_baoab!(psB, f, dt; γ=γ, T=Tt, kB=1.0, rng=rngB)
        end
        @test psA.positions == psB.positions
        @test psA.velocities == psB.velocities
    end
end
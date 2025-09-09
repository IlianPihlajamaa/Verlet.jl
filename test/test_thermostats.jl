using Test
using Random
using Verlet
using StaticArrays

@testset "thermostats.jl" begin
    @testset "instantaneous_temperature & degrees_of_freedom" begin
    N, D = 10, 3
    pos = [SVector{D,Float64}(0.0,0.0,0.0) for _ in 1:N]
    vmag = sqrt(2 * 0.5 * N * D * 1.0 * 2.0 / sum(ones(N)))
    vel = [SVector{D,Float64}(vmag,0.0,0.0) for _ in 1:N]
    ps = ParticleSystem{D,Float64}(pos, vel, ones(N))
    Tt, kB = 2.0, 1.0
    dof = N*D
    @test degrees_of_freedom(ps) == dof
    Ti = instantaneous_temperature(ps; kB=kB)
    @test isapprox(Ti, Tt; rtol=1e-12, atol=0)
    end

    @testset "velocity_rescale! matches target T" begin
    N, D = 8, 3
    pos = [SVector{D,Float64}(0.0,0.0,0.0) for _ in 1:N]
    vel = [SVector{D,Float64}(randn(3)...) for _ in 1:N]
    ps = ParticleSystem{D,Float64}(pos, vel, ones(N))
    Tt = 3.0
    velocity_rescale!(ps, Tt; kB=1.0)
    @test isapprox(instantaneous_temperature(ps; kB=1.0), Tt; rtol=1e-12, atol=0)
    end

    @testset "BAOAB: γ=0 reduces to VV trajectory (single step)" begin
        # Compare one BAOAB step with γ=0 against velocity_verlet! for a harmonic oscillator.
        using LinearAlgebra, StaticArrays
    N, D = 4, 3
    k = 1.0
    ho_forces(R) = map(x -> -k * x, R)
    pos = [SVector{D,Float64}(randn(3)...) for _ in 1:N]
    vel = [SVector{D,Float64}(randn(3)...) for _ in 1:N]
    ps1 = ParticleSystem{D,Float64}(copy(pos), copy(vel), ones(N))
    ps2 = ParticleSystem{D,Float64}(copy(pos), copy(vel), ones(N))
    dt = 0.01
    rng = MersenneTwister(42)
    langevin_baoab!(ps1, ho_forces, dt; γ=0.0, T=0.0, kB=1.0, rng=rng)
    velocity_verlet!(ps2, ho_forces, dt)
    @test all(isapprox.(ps1.positions, ps2.positions; rtol=1e-12, atol=1e-12))
    @test all(isapprox.(ps1.velocities, ps2.velocities; rtol=1e-12, atol=1e-12))
    end

    @testset "BAOAB: temperature control converges to target (statistical)" begin
        # Stochastic test with tolerance: run several steps and check mean T near target.
        N, D = 5000, 3
        pos = [SVector{D,Float64}(0.0,0.0,0.0) for _ in 1:N]
        vel = [SVector{D,Float64}(randn(3)...) for _ in 1:N]
        ps = ParticleSystem{D,Float64}(pos, vel, ones(N))
        Tt, γ, dt = 1.5, 1.0, 0.005
        rng = MersenneTwister(123)
        steps = 2000
        for _ in 1:steps
            ho_forces(R) = -0.1 * R
            langevin_baoab!(ps, ho_forces, dt; γ=γ, T=Tt, kB=1.0, rng=rng)
        end
        Tinst = instantaneous_temperature(ps; kB=1.0)
        @test isapprox(Tinst, Tt; rtol=0.05)
    end

    @testset "BAOAB reproducibility with fixed RNG" begin
        N, D = 16, 3
        pos = [SVector{D,Float64}(randn(3)...) for _ in 1:N]
        vel = [SVector{D,Float64}(randn(3)...) for _ in 1:N]
        psA = ParticleSystem{D,Float64}(copy(pos), copy(vel), ones(N))
        psB = ParticleSystem{D,Float64}(copy(pos), copy(vel), ones(N))
        Tt, γ, dt = 2.0, 0.7, 0.01
        rngA = MersenneTwister(2025)
        rngB = MersenneTwister(2025)
        steps = 50
        for _ in 1:steps
            f(R) = map(x -> -0.2 * x, R)
            langevin_baoab!(psA, f, dt; γ=γ, T=Tt, kB=1.0, rng=rngA)
            langevin_baoab!(psB, f, dt; γ=γ, T=Tt, kB=1.0, rng=rngB)
        end
        @test psA.positions == psB.positions
        @test psA.velocities == psB.velocities
    end
end
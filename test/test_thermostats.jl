using Test
using Random
using Verlet
using StaticArrays

struct ThermoHarmonic end

function Verlet.Core.compute_forces!(::ThermoHarmonic, sys::System)
    @inbounds for i in 1:natoms(sys)
        sys.forces[i] += -sys.positions[i]
    end
    return sys
end

@testset "thermostats.jl" begin
    @testset "instantaneous_temperature & degrees_of_freedom" begin
        N, D = 10, 3
        positions = [SVector{D,Float64}(0.0,0.0,0.0) for _ in 1:N]
        vmag = sqrt(2 * 0.5 * N * D * 1.0 * 2.0 / sum(ones(N)))
        velocities = [SVector{D,Float64}(vmag,0.0,0.0) for _ in 1:N]
        forces = [SVector{D,Float64}(0.0,0.0,0.0) for _ in 1:N]
        masses = ones(N)
        box = CubicBox(10.0)
        types = ones(Int, N)
        type_names = Dict(1 => :A)
        sys = System(positions, velocities, forces, masses, box, types, type_names)
        Tt, kB = 2.0, 1.0
        dof = N*D
        @test degrees_of_freedom(sys) == dof
        Ti = instantaneous_temperature(sys; kB=kB)
        @test isapprox(Ti, Tt; rtol=1e-12, atol=0)
    end

    @testset "velocity_rescale! matches target T" begin
        N, D = 8, 3
        positions = [SVector{D,Float64}(0.0,0.0,0.0) for _ in 1:N]
        velocities = [SVector{D,Float64}(randn(3)...) for _ in 1:N]
        forces = [SVector{D,Float64}(0.0,0.0,0.0) for _ in 1:N]
        masses = ones(N)
        box = CubicBox(10.0)
        types = ones(Int, N)
        type_names = Dict(1 => :A)
        sys = System(positions, velocities, forces, masses, box, types, type_names)
        Tt = 3.0
        velocity_rescale!(sys, Tt; kB=1.0)
        @test isapprox(instantaneous_temperature(sys; kB=1.0), Tt; rtol=1e-12, atol=0)
    end

    @testset "BAOAB: γ=0 reduces to VV trajectory (single step)" begin
        # Compare one BAOAB step with γ=0 against VelocityVerlet for a harmonic oscillator.
        using LinearAlgebra, StaticArrays
        N, D = 4, 3
        k = 1.0
        ho_forces(R) = map(x -> -k * x, R)
        positions = [SVector{D,Float64}(randn(3)...) for _ in 1:N]
        velocities = [SVector{D,Float64}(randn(3)...) for _ in 1:N]
        forces = [SVector{D,Float64}(0.0,0.0,0.0) for _ in 1:N]
        masses = ones(N)
        box = CubicBox(10.0)
        types = ones(Int, N)
        type_names = Dict(1 => :A)
        sys1 = System(copy(positions), copy(velocities), copy(forces), masses, box, types, type_names)
        sys2 = System(copy(positions), copy(velocities), copy(forces), masses, box, types, type_names;
                      forcefield=ForceField((ThermoHarmonic(),)))
        dt = 0.01
        rng = MersenneTwister(42)
        langevin_baoab!(sys1, ho_forces, dt; γ=0.0, temp=0.0, kB=1.0, rng=rng)
        integrate!(VelocityVerlet(dt), sys2, 1)
        @test all(isapprox.(sys1.positions, sys2.positions; rtol=1e-12, atol=1e-12))
        @test all(isapprox.(sys1.velocities, sys2.velocities; rtol=1e-12, atol=1e-12))
    end

    @testset "BAOAB: temperature control converges to target (statistical)" begin
        # Stochastic test with tolerance: run several steps and check mean T near target.
        N, D = 5000, 3
        positions = [SVector{D,Float64}(0.0,0.0,0.0) for _ in 1:N]
        velocities = [SVector{D,Float64}(randn(3)...) for _ in 1:N]
        forces = [SVector{D,Float64}(0.0,0.0,0.0) for _ in 1:N]
        masses = ones(N)
        box = CubicBox(10.0)
        types = ones(Int, N)
        type_names = Dict(1 => :A)
        sys = System(positions, velocities, forces, masses, box, types, type_names)
        Tt, γ, dt = 1.5, 1.0, 0.005
        rng = MersenneTwister(123)
        steps = 2000
        for _ in 1:steps
            ho_forces(R) = -0.1 * R
            langevin_baoab!(sys, ho_forces, dt; γ=γ, temp=Tt, kB=1.0, rng=rng)
        end
        Tinst = instantaneous_temperature(sys; kB=1.0)
        @test isapprox(Tinst, Tt; rtol=0.05)
    end

    @testset "BAOAB reproducibility with fixed RNG" begin
        N, D = 16, 3
        positions = [SVector{D,Float64}(randn(3)...) for _ in 1:N]
        velocities = [SVector{D,Float64}(randn(3)...) for _ in 1:N]
        forces = [SVector{D,Float64}(0.0,0.0,0.0) for _ in 1:N]
        masses = ones(N)
        box = CubicBox(10.0)
        types = ones(Int, N)
        type_names = Dict(1 => :A)
        sysA = System(copy(positions), copy(velocities), copy(forces), masses, box, types, type_names)
        sysB = System(copy(positions), copy(velocities), copy(forces), masses, box, types, type_names)
        Tt, γ, dt = 2.0, 0.7, 0.01
        rngA = MersenneTwister(2025)
        rngB = MersenneTwister(2025)
        steps = 50
        for _ in 1:steps
            f(R) = map(x -> -0.2 * x, R)
            langevin_baoab!(sysA, f, dt; γ=γ, temp=Tt, kB=1.0, rng=rngA)
            langevin_baoab!(sysB, f, dt; γ=γ, temp=Tt, kB=1.0, rng=rngB)
        end
        @test sysA.positions == sysB.positions
        @test sysA.velocities == sysB.velocities
    end
end

using Test
using Random
using Verlet
using StaticArrays

struct ThermoHarmonic
    k::Float64
end

function Verlet.Core.compute_forces!(pot::ThermoHarmonic, sys::System)
    @inbounds for i in 1:natoms(sys)
        sys.forces[i] += -pot.k * sys.positions[i]
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
        N, D = 4, 3
        positions = [SVector{D,Float64}(randn(3)...) for _ in 1:N]
        velocities = [SVector{D,Float64}(randn(3)...) for _ in 1:N]
        forces = [SVector{D,Float64}(0.0,0.0,0.0) for _ in 1:N]
        masses = ones(N)
        box = CubicBox(10.0)
        types = ones(Int, N)
        type_names = Dict(1 => :A)
        ff = ForceField((ThermoHarmonic(1.0),))
        sys1 = System(copy(positions), copy(velocities), copy(forces), masses, box, types, type_names; forcefield=ff)
        sys2 = System(copy(positions), copy(velocities), copy(forces), masses, box, types, type_names; forcefield=ff)
        dt = 0.01
        rng = MersenneTwister(42)
        integrate!(LangevinBAOAB(dt; γ=0.0, temp=0.0, rng=rng), sys1, 1)
        integrate!(VelocityVerlet(dt), sys2, 1)
        @test all(isapprox.(sys1.positions, sys2.positions; rtol=1e-12, atol=1e-12))
        @test all(isapprox.(sys1.velocities, sys2.velocities; rtol=1e-12, atol=1e-12))
    end

    @testset "BAOAB: temperature control converges to target (statistical)" begin
        N, D = 5000, 3
        positions = [SVector{D,Float64}(0.0,0.0,0.0) for _ in 1:N]
        velocities = [SVector{D,Float64}(randn(3)...) for _ in 1:N]
        forces = [SVector{D,Float64}(0.0,0.0,0.0) for _ in 1:N]
        masses = ones(N)
        box = CubicBox(10.0)
        types = ones(Int, N)
        type_names = Dict(1 => :A)
        ff = ForceField((ThermoHarmonic(0.1),))
        sys = System(positions, velocities, forces, masses, box, types, type_names; forcefield=ff)
        Tt, γ, dt = 1.5, 1.0, 0.005
        rng = MersenneTwister(123)
        steps = 2000
        integrate!(LangevinBAOAB(dt; γ=γ, temp=Tt, rng=rng), sys, steps)
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
        ff = ForceField((ThermoHarmonic(0.2),))
        sysA = System(copy(positions), copy(velocities), copy(forces), masses, box, types, type_names; forcefield=ff)
        sysB = System(copy(positions), copy(velocities), copy(forces), masses, box, types, type_names; forcefield=ff)
        Tt, γ, dt = 2.0, 0.7, 0.01
        rngA = MersenneTwister(2025)
        rngB = MersenneTwister(2025)
        steps = 50
        integrate!(LangevinBAOAB(dt; γ=γ, temp=Tt, rng=rngA), sysA, steps)
        integrate!(LangevinBAOAB(dt; γ=γ, temp=Tt, rng=rngB), sysB, steps)
        @test sysA.positions == sysB.positions
        @test sysA.velocities == sysB.velocities
    end
end

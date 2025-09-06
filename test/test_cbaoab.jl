using Test
using Random
using LinearAlgebra, StaticArrays
using Verlet

@testset "cBAOAB integrator with constraints" begin
@testset "cBAOAB preserves constraints (zero forces)" begin
    N, D = 3, 3
    r0 = 1.0
    positions = [@SVector zeros(D) for _ in 1:N]
    positions[2] = @SVector [r0, 0.0, 0.0]
    velocities = [@SVector zeros(D) for _ in 1:N]
    masses = ones(N)
    ps = ParticleSystem(positions, velocities, masses)
    cons = DistanceConstraints([(1,2)], [r0]; tol=1e-10, maxiter=200)

    forces(R) = [@SVector zeros(D) for _ in R]
    dt, γ, T = 0.01, 1.0, 0.5
    rng = MersenneTwister(1234)

    # Warmup and integrate
    for _ in 1:1000
        langevin_baoab_constrained!(ps, forces, dt, cons; γ=γ, T=T, rng=rng)
    end

    d = ps.positions[1] - ps.positions[2]
    @test isapprox(norm(d), r0; atol=1e-6)

    # Velocity residual should be tiny
    (; maxC, maxCd) = constraint_residuals(ps, cons)
    @test maxC ≤ 1e-8
    @test maxCd ≤ 1e-8
end

@testset "cBAOAB temperature sane and finite" begin
    N, D = 8, 3
    positions = [@SVector randn(D) for _ in 1:N]
    velocities = [@SVector zeros(D) for _ in 1:N]
    masses = ones(N)
    ps = ParticleSystem(positions, velocities, masses)
    cons = DistanceConstraints([(1,2)], [1.0])
    forces(R) = [@SVector zeros(D) for _ in R]
    dt, γ, T = 0.005, 2.0, 1.2
    rng = MersenneTwister(7)

    # Run; measure running average temperature
    accT = 0.0
    steps = 5000
    for n in 1:steps
        langevin_baoab_constrained!(ps, forces, dt, cons; γ=γ, T=T, rng=rng)
        accT += instantaneous_temperature(ps; kB=1.0)
    end
    Tbar = accT / steps
    @test isfinite(Tbar)
    # Loose bound (small system; not asserting equality to T)
    @test 0.2 ≤ Tbar ≤ 3.0
end

@testset "constraint_residuals reports zero at exact satisfaction" begin
    positions = [@SVector([0.0, 0.0, 0.0]), @SVector([1.0, 0.0, 0.0])]
    velocities = [@SVector(zeros(3)) for _ in 1:2]
    masses = ones(2)
    ps = ParticleSystem(positions, velocities, masses)
    cons = DistanceConstraints([(1,2)], [1.0])
    (; maxC, rmsC, maxCd, rmsCd) = constraint_residuals(ps, cons)
    @test maxC == 0.0
    @test rmsC == 0.0
    @test maxCd == 0.0
    @test rmsCd == 0.0
end
end
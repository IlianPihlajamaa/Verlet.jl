
using Test
using Verlet
using StaticArrays


@testset "Constraints" begin
    @testset "DistanceConstraints basics" begin
    N, D = 2, 3
    r0 = 1.25
    pos = [SVector{D,Float64}(0.0,0.0,0.0), SVector{D,Float64}(r0+0.1,0.0,0.0)]
    vel = [SVector{D,Float64}(0.0,0.0,0.0) for _ in 1:N]
    ps = ParticleSystem{D,Float64}(pos, vel, ones(N))
    cons = DistanceConstraints([(1,2)], [r0]; tol=1e-10, maxiter=100)
    apply_shake!(ps, cons, 0.01)
    d = ps.positions[1] - ps.positions[2]
    @test isapprox(norm(d), r0; rtol=0, atol=1e-8)
    end

    @testset "RATTLE projects velocities" begin
    N, D = 2, 3
    r0 = 1.0
    pos = [SVector{D,Float64}(0.0,0.0,0.0), SVector{D,Float64}(r0,0.0,0.0)]
    vel = [SVector{D,Float64}(+0.3,0.0,0.0), SVector{D,Float64}(-0.1,0.0,0.0)]
    ps = ParticleSystem{D,Float64}(pos, vel, ones(N))
    cons = DistanceConstraints([(1,2)], [r0])
    apply_rattle!(ps, cons)
    d = ps.positions[1] - ps.positions[2]
    vrel = ps.velocities[1] - ps.velocities[2]
    @test isapprox(dot(d, vrel), 0.0; atol=1e-10)
    end

    @testset "Constrained VV holds bond length" begin
        N, D = 2, 3
        r0 = 0.75
        pos = [SVector{D,Float64}(0.0,0.0,0.0), SVector{D,Float64}(r0,0.0,0.0)]
        vel = [SVector{D,Float64}(0.01,0.01,0.01) for _ in 1:N]
        ps = ParticleSystem{D,Float64}(pos, vel, ones(N))
        cons = DistanceConstraints([(1,2)], [r0]; tol=1e-10, maxiter=100)
        dt = 0.02
        forces(R) = [SVector{D,Float64}(0.0,0.0,0.0) for _ in 1:N]
        for _ in 1:500
            velocity_verlet_shake_rattle!(ps, forces, dt, cons)
        end
        d = ps.positions[1] - ps.positions[2]
        @test isapprox(norm(d), r0; atol=1e-7)
    end

    @testset "degrees_of_freedom with constraints and COM" begin
    N, D = 5, 3
    pairs = [(1,2), (3,4)]
    r0s = [1.0, 1.5]
    cons = DistanceConstraints(pairs, r0s)
    pos = [SVector{D,Float64}(0.0,0.0,0.0) for _ in 1:N]
    vel = [SVector{D,Float64}(0.0,0.0,0.0) for _ in 1:N]
    ps = ParticleSystem{D,Float64}(pos, vel, ones(N))
    @test degrees_of_freedom(ps; constraints=cons, remove_com=false) == N*D - length(pairs)
    @test degrees_of_freedom(ps; constraints=cons, remove_com=true) == N*D - length(pairs) - D
    end

    @testset "remove_com_motion! removes mass-weighted COM velocity" begin
    N, D = 3, 3
    pos = [SVector{D,Float64}(0.0,0.0,0.0) for _ in 1:N]
    vel = [SVector{D,Float64}(1.0,1.0,1.0) for _ in 1:N]
    ps = ParticleSystem{D,Float64}(pos, vel, [1.0, 2.0, 3.0])
    remove_com_motion!(ps; which=:velocity)

    Vcom = sum(ps.masses[i] * ps.velocities[i] for i in 1:N) / sum(ps.masses)

    @test all(isapprox.(Vcom, 0.0; atol=1e-14))
    end
end
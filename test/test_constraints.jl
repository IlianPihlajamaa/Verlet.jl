using Test
using Verlet


@testset "Constraints" begin
    @testset "DistanceConstraints basics" begin
        N, D = 2, 3
        r0 = 1.25
        ps = ParticleSystem(zeros(N,D), zeros(N,D), ones(N))
        ps.positions[1,1] = 0.0
        ps.positions[2,1] = r0 + 0.1
        cons = DistanceConstraints([(1,2)], [r0]; tol=1e-10, maxiter=100)
        apply_shake!(ps, cons, 0.01)
        d = ps.positions[1,:] .- ps.positions[2,:]
        @test isapprox(norm(d), r0; rtol=0, atol=1e-8)
    end

    @testset "RATTLE projects velocities" begin
        N, D = 2, 3
        r0 = 1.0
        ps = ParticleSystem(zeros(N,D), zeros(N,D), ones(N))
        ps.positions[1,1] = 0.0
        ps.positions[2,1] = r0
        cons = DistanceConstraints([(1,2)], [r0])
        ps.velocities[1,1] = +0.3
        ps.velocities[2,1] = -0.1
        apply_rattle!(ps, cons)
        d = ps.positions[1,:] .- ps.positions[2,:]
        vrel = ps.velocities[1,:] .- ps.velocities[2,:]
        @test isapprox(dot(d, vrel), 0.0; atol=1e-10)
    end

    @testset "Constrained VV holds bond length" begin
        N, D = 2, 3
        r0 = 0.75
        ps = ParticleSystem(zeros(N,D), zeros(N,D), ones(N))
        ps.positions[2,1] = r0
        ps.velocities .= 0.01
        cons = DistanceConstraints([(1,2)], [r0]; tol=1e-10, maxiter=100)
        dt = 0.02
        forces(R) = zero(R)
        for _ in 1:500
            velocity_verlet_shake_rattle!(ps, forces, dt, cons)
        end
        d = ps.positions[1,:] .- ps.positions[2,:]
        @test isapprox(norm(d), r0; atol=1e-7)
    end

    @testset "degrees_of_freedom with constraints and COM" begin
        N, D = 5, 3
        pairs = [(1,2), (3,4)]
        r0s = [1.0, 1.5]
        cons = DistanceConstraints(pairs, r0s)
        ps = ParticleSystem(zeros(N,D), zeros(N,D), ones(N))
        @test degrees_of_freedom(ps; constraints=cons, remove_com=false) == N*D - length(pairs)
        @test degrees_of_freedom(ps; constraints=cons, remove_com=true) == N*D - length(pairs) - D
    end

    @testset "remove_com_motion! removes mass-weighted COM velocity" begin
        N, D = 3, 3
        ps = ParticleSystem(zeros(N,D), zeros(N,D), [1.0, 2.0, 3.0])
        ps.velocities .= 1.0
        remove_com_motion!(ps; which=:velocity)
        Vcom = (sum(ps.masses .* ps.velocities[:,1])) / sum(ps.masses)
        @test isapprox(Vcom, 0.0; atol=1e-14)
    end
end
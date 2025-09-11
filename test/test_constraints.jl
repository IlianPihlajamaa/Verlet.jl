using Test
using Verlet
using StaticArrays


@testset "Constraints" begin
    @testset "DistanceConstraints basics" begin
        N, D = 2, 3
        r0 = 1.25
        positions = [SVector{D,Float64}(0.0,0.0,0.0), SVector{D,Float64}(r0+0.1,0.0,0.0)]
        velocities = [SVector{D,Float64}(0.0,0.0,0.0) for _ in 1:N]
        forces = [SVector{D,Float64}(0.0,0.0,0.0) for _ in 1:N]
        masses = ones(N)
        box = CubicBox(10.0)
        types = ones(Int, N)
        type_names = Dict(1 => :A)
        sys = System(positions, velocities, forces, masses, box, types, type_names)
        cons = DistanceConstraints([(1,2)], [r0]; tol=1e-10, maxiter=100)
        apply_shake!(sys, cons, 0.01)
        d = sys.positions[1] - sys.positions[2]
        @test isapprox(norm(d), r0; rtol=0, atol=1e-8)
    end

    @testset "RATTLE projects velocities" begin
        N, D = 2, 3
        r0 = 1.0
        positions = [SVector{D,Float64}(0.0,0.0,0.0), SVector{D,Float64}(r0,0.0,0.0)]
        velocities = [SVector{D,Float64}(+0.3,0.0,0.0), SVector{D,Float64}(-0.1,0.0,0.0)]
        forces = [SVector{D,Float64}(0.0,0.0,0.0) for _ in 1:N]
        masses = ones(N)
        box = CubicBox(10.0)
        types = ones(Int, N)
        type_names = Dict(1 => :A)
        sys = System(positions, velocities, forces, masses, box, types, type_names)
        cons = DistanceConstraints([(1,2)], [r0])
        apply_rattle!(sys, cons)
        d = sys.positions[1] - sys.positions[2]
        vrel = sys.velocities[1] - sys.velocities[2]
        @test isapprox(dot(d, vrel), 0.0; atol=1e-10)
    end

    @testset "Constrained VV holds bond length" begin
        N, D = 2, 3
        r0 = 0.75
        positions = [SVector{D,Float64}(0.0,0.0,0.0), SVector{D,Float64}(r0,0.0,0.0)]
        velocities = [SVector{D,Float64}(0.01,0.01,0.01) for _ in 1:N]
        forces_storage = [SVector{D,Float64}(0.0,0.0,0.0) for _ in 1:N]
        masses = ones(N)
        box = CubicBox(10.0)
        types = ones(Int, N)
        type_names = Dict(1 => :A)
        sys = System(positions, velocities, forces_storage, masses, box, types, type_names)
        cons = DistanceConstraints([(1,2)], [r0]; tol=1e-10, maxiter=100)
        dt = 0.02
        forces_func(R) = [SVector{D,Float64}(0.0,0.0,0.0) for _ in 1:N]
        for _ in 1:500
            velocity_verlet_shake_rattle!(sys, forces_func, dt, cons)
        end
        d = sys.positions[1] - sys.positions[2]
        @test isapprox(norm(d), r0; atol=1e-7)
    end

    @testset "degrees_of_freedom with constraints and COM" begin
        N, D = 5, 3
        pairs = [(1,2), (3,4)]
        r0s = [1.0, 1.5]
        cons = DistanceConstraints(pairs, r0s)
        positions = [SVector{D,Float64}(0.0,0.0,0.0) for _ in 1:N]
        velocities = [SVector{D,Float64}(0.0,0.0,0.0) for _ in 1:N]
        forces = [SVector{D,Float64}(0.0,0.0,0.0) for _ in 1:N]
        masses = ones(N)
        box = CubicBox(10.0)
        types = ones(Int, N)
        type_names = Dict(1 => :A)
        sys = System(positions, velocities, forces, masses, box, types, type_names)
        @test degrees_of_freedom(sys; constraints=cons, remove_com=false) == N*D - length(pairs)
        @test degrees_of_freedom(sys; constraints=cons, remove_com=true) == N*D - length(pairs) - D
    end

    @testset "remove_com_motion! removes mass-weighted COM velocity" begin
        N, D = 3, 3
        positions = [SVector{D,Float64}(0.0,0.0,0.0) for _ in 1:N]
        velocities = [SVector{D,Float64}(1.0,1.0,1.0) for _ in 1:N]
        forces = [SVector{D,Float64}(0.0,0.0,0.0) for _ in 1:N]
        masses = [1.0, 2.0, 3.0]
        box = CubicBox(10.0)
        types = ones(Int, N)
        type_names = Dict(1 => :A)
        sys = System(positions, velocities, forces, masses, box, types, type_names)
        remove_com_motion!(sys; which=:velocity)

        Vcom = sum(sys.masses[i] * sys.velocities[i] for i in 1:N) / sum(sys.masses)

        @test all(isapprox.(Vcom, 0.0; atol=1e-14))
    end
end
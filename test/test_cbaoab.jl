using Test
using Random
using LinearAlgebra, StaticArrays
using Verlet

struct ZeroForce end

function Verlet.Core.compute_forces!(::ZeroForce, sys::System)
    fill!(sys.forces, zero(sys.forces[1]))
    return sys
end

@testset "cBAOAB integrator with constraints" begin
    @testset "cBAOAB preserves constraints (zero forces)" begin
        N, D = 3, 3
        r0 = 1.0
        positions = [@SVector zeros(D) for _ in 1:N]
        positions[2] = @SVector [r0, 0.0, 0.0]
        velocities = [@SVector zeros(D) for _ in 1:N]
        forces = [@SVector zeros(D) for _ in 1:N]
        masses = ones(N)
        box = CubicBox(10.0)
        types = ones(Int, N)
        type_names = Dict(1 => :A)
        sys = System(positions, velocities, forces, masses, box, types, type_names;
                     forcefield=ForceField((ZeroForce(),)))
        cons = DistanceConstraints([(1,2)], [r0]; tol=1e-10, maxiter=200)

        dt, γ, T = 0.01, 1.0, 0.5
        rng = MersenneTwister(1234)
        integrator = LangevinBAOABConstrained(dt, cons; γ=γ, temp=T, rng=rng)

        integrate!(integrator, sys, 1000)

        d = sys.positions[1] - sys.positions[2]
        @test isapprox(norm(d), r0; atol=1e-6)

        (; maxC, maxCd) = constraint_residuals(sys, cons)
        @test maxC ≤ 1e-8
        @test maxCd ≤ 1e-8
    end

    @testset "cBAOAB temperature sane and finite" begin
        N, D = 8, 3
        positions = [@SVector randn(D) for _ in 1:N]
        velocities = [@SVector zeros(D) for _ in 1:N]
        forces = [@SVector zeros(D) for _ in 1:N]
        masses = ones(N)
        box = CubicBox(10.0)
        types = ones(Int, N)
        type_names = Dict(1 => :A)
        sys = System(positions, velocities, forces, masses, box, types, type_names;
                     forcefield=ForceField((ZeroForce(),)))
        cons = DistanceConstraints([(1,2)], [1.0])
        dt, γ, T = 0.005, 2.0, 1.2
        rng = MersenneTwister(7)

        integrator = LangevinBAOABConstrained(dt, cons; γ=γ, temp=T, rng=rng)
        accT = 0.0
        steps = 5000
        integrate!(integrator, sys, steps;
                   callback = (sys, step, _) -> begin
                       accT += instantaneous_temperature(sys; kB=1.0)
                       return nothing
                   end)
        Tbar = accT / steps
        @test isfinite(Tbar)
        @test 0.2 ≤ Tbar ≤ 3.0
    end

    @testset "constraint_residuals reports zero at exact satisfaction" begin
        positions = [@SVector([0.0, 0.0, 0.0]), @SVector([1.0, 0.0, 0.0])]
        velocities = [@SVector(zeros(3)) for _ in 1:2]
        forces = [@SVector(zeros(3)) for _ in 1:2]
        masses = ones(2)
        box = CubicBox(10.0)
        types = ones(Int, 2)
        type_names = Dict(1 => :A)
        sys = System(positions, velocities, forces, masses, box, types, type_names;
                     forcefield=ForceField((ZeroForce(),)))
        cons = DistanceConstraints([(1,2)], [1.0])
        (; maxC, rmsC, maxCd, rmsCd) = constraint_residuals(sys, cons)
        @test maxC == 0.0
        @test rmsC == 0.0
        @test maxCd == 0.0
        @test rmsCd == 0.0
    end

    @testset "integrator matches manual constrained BAOAB" begin
        D = 3
        positions = [
            @SVector([0.0, 0.0, 0.0]),
            @SVector([1.0, 0.0, 0.0]),
            @SVector([0.0, 1.0, 0.0]),
            @SVector([0.0, 0.0, 1.0]),
        ]
        N = length(positions)
        velocities = [@SVector zeros(D) for _ in 1:N]
        forces = [@SVector zeros(D) for _ in 1:N]
        masses = ones(N)
        box = CubicBox(8.0)
        types = ones(Int, N)
        type_names = Dict(1 => :A)
        ff = ForceField((ZeroForce(),))
        sys1 = System(deepcopy(positions), deepcopy(velocities), deepcopy(forces), deepcopy(masses), box, types, type_names; forcefield=ff)
        sys2 = deepcopy(sys1)
        cons1 = DistanceConstraints([(1,2)], [1.0])
        cons2 = DistanceConstraints([(1,2)], [1.0])
        dt, γ, T = 0.005, 1.5, 0.8
        rng1 = MersenneTwister(2024)
        rng2 = MersenneTwister(2024)
        steps = 7

        integrate!(LangevinBAOABConstrained(dt, cons1; γ=γ, temp=T, rng=rng1), sys1, steps)
        integrate!(LangevinBAOABConstrained(dt, cons2; γ=γ, temp=T, rng=rng2), sys2, steps)

        pos_err = maximum(norm.(sys1.positions .- sys2.positions))
        vel_err = maximum(norm.(sys1.velocities .- sys2.velocities))
        @test pos_err ≤ 1e-10
        @test vel_err ≤ 1e-10

        res = constraint_residuals(sys1, cons1)
        @test res.maxC ≤ 1e-8
        @test res.maxCd ≤ 1e-8
    end
end

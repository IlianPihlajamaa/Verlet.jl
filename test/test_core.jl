using Test
using Verlet
using StaticArrays
using LinearAlgebra, StaticArrays

const Dims = 3

@testset "Core functionality" begin
    @testset "System construction" begin
        positions = [SVector{Dims,Float64}(0.0,0.0,0.0)]
        velocities = [SVector{Dims,Float64}(1.0,0.0,0.0)]
        forces = [SVector{Dims,Float64}(0.0,0.0,0.0)]
        masses = [1.0]
        box = CubicBox(10.0)
        types = [1]
        type_names = Dict(1 => :A)
        sys = System(positions, velocities, forces, masses, box, types, type_names)
        @test length(sys.positions) == 1
        @test length(sys.velocities) == 1
        @test sys.masses == [1.0]
    end
    @testset "Velocity Verlet" begin
        @testset "Velocity Verlet free particle" begin
            forces_func(r) = fill(SVector{Dims,Float64}(0.0,0.0,0.0), length(r))
            positions = [SVector{Dims,Float64}(0.0,0.0,0.0)]
            velocities = [SVector{Dims,Float64}(1.0,0.0,0.0)]
            forces = [SVector{Dims,Float64}(0.0,0.0,0.0)]
            masses = [1.0]
            box = CubicBox(10.0)
            types = [1]
            type_names = Dict(1 => :A)
            sys = System(positions, velocities, forces, masses, box, types, type_names)
            velocity_verlet!(sys, forces_func, 0.1)
            @test isapprox(sys.positions[1][1], 0.1; atol=1e-12)
            @test isapprox(sys.velocities[1][1], 1.0; atol=1e-12)
        end

        @testset "Velocity Verlet harmonic oscillator" begin
            forces_func(r) = map(x -> -x, r)
            positions = [SVector{Dims,Float64}(1.0,0.0,0.0)]
            velocities = [SVector{Dims,Float64}(0.0,0.0,0.0)]
            forces = [SVector{Dims,Float64}(0.0,0.0,0.0)]
            masses = [1.0]
            box = CubicBox(10.0)
            types = [1]
            type_names = Dict(1 => :A)
            sys = System(positions, velocities, forces, masses, box, types, type_names)
            velocity_verlet!(sys, forces_func, 0.1)
            @test sys.positions[1][1] < 1.0   # should move left
        end
    end

    @testset "integrate! driver" begin
        forces_func(r) = map(x -> -x, r)
        function make_system(; pos=SVector{Dims,Float64}(1.0, 0.0, 0.0),
                               vel=SVector{Dims,Float64}(0.0, 1.0, 0.0))
            positions = [pos]
            velocities = [vel]
            forces = [SVector{Dims,Float64}(0.0, 0.0, 0.0)]
            masses = [1.0]
            box = CubicBox(10.0)
            types = [1]
            type_names = Dict(1 => :A)
            return System(positions, velocities, forces, masses, box, types, type_names)
        end
        dt = 0.1

        @testset "matches direct loop" begin
            sys_integrated = make_system()
            sys_reference = deepcopy(sys_integrated)
            integrate!(velocity_verlet!, sys_integrated, forces_func, dt, 5)
            for _ in 1:5
                velocity_verlet!(sys_reference, forces_func, dt)
            end
            @test all(isapprox.(sys_integrated.positions, sys_reference.positions; atol=1e-12))
            @test all(isapprox.(sys_integrated.velocities, sys_reference.velocities; atol=1e-12))
        end

        @testset "zero steps no-op" begin
            sys_zero = make_system()
            sys_copy = deepcopy(sys_zero)
            calls = Ref(0)
            integrate!(velocity_verlet!, sys_zero, forces_func, dt, 0; callback = (sys, step) -> (calls[] += 1))
            @test sys_zero.positions == sys_copy.positions
            @test sys_zero.velocities == sys_copy.velocities
            @test calls[] == 0
        end

        @testset "callback early stop" begin
            sys_cb = make_system()
            sys_reference = deepcopy(sys_cb)
            last_step = Ref(0)
            integrate!(velocity_verlet!, sys_cb, forces_func, dt, 10;
                       callback = (sys, step) -> begin
                           last_step[] = step
                           return step == 3 ? false : nothing
                       end)
            for _ in 1:3
                velocity_verlet!(sys_reference, forces_func, dt)
            end
            @test last_step[] == 3
            @test all(isapprox.(sys_cb.positions, sys_reference.positions; atol=1e-12))
            @test all(isapprox.(sys_cb.velocities, sys_reference.velocities; atol=1e-12))
        end

        function drift_with_offset!(sys::System, forces, dt, offset; scale::Real=1.0)
            forces(sys.positions)
            for i in eachindex(sys.positions)
                sys.positions[i] += offset * (dt * scale)
            end
            return sys
        end
        zero_forces(R) = fill(SVector{Dims,Float64}(0.0, 0.0, 0.0), length(R))

        @testset "forwards extra args" begin
            positions = [SVector{Dims,Float64}(0.0, 0.0, 0.0)]
            velocities = [SVector{Dims,Float64}(0.0, 0.0, 0.0)]
            forces = [SVector{Dims,Float64}(0.0, 0.0, 0.0)]
            masses = [1.0]
            box = CubicBox(5.0)
            types = [1]
            type_names = Dict(1 => :A)
            sys_extra = System(positions, velocities, forces, masses, box, types, type_names)
            offset = SVector{Dims,Float64}(1.0, 0.0, 0.0)
            integrate!(drift_with_offset!, sys_extra, zero_forces, 0.5, 4, offset; scale=2.0)
            @test sys_extra.positions[1] == SVector{Dims,Float64}(4.0, 0.0, 0.0)
            @test sys_extra.velocities[1] == SVector{Dims,Float64}(0.0, 0.0, 0.0)
        end

        @testset "rejects negative steps" begin
            sys_neg = make_system()
            @test_throws ArgumentError integrate!(velocity_verlet!, sys_neg, forces_func, dt, -2)
        end
    end

    @testset "Velocity Verlet energy conservation" begin
        function ho_forces(R; return_potential::Bool=false)
            F = [-r for r in R]
            if return_potential
                U = 0.5 * sum(norm(r)^2 for r in R)
                return F, U
            end
            return F
        end

        positions = [SVector{Dims,Float64}(1.0, 0.0, 0.0)]
        velocities = [SVector{Dims,Float64}(0.0, 1.0, 0.0)]
        forces = [SVector{Dims,Float64}(0.0, 0.0, 0.0)]
        masses = [1.0]
        box = CubicBox(10.0)
        types = [1]
        type_names = Dict(1 => :A)
        sys = System(positions, velocities, forces, masses, box, types, type_names)

        (_, U0) = ho_forces(sys.positions; return_potential=true)
        KE0 = kinetic_energy(sys)
        E0 = KE0 + U0

        energies = Float64[E0]
        dt = 0.002
        steps = 2000
        integrate!(velocity_verlet!, sys, ho_forces, dt, steps;
                   callback = (sys, step) -> begin
                       (_, U) = ho_forces(sys.positions; return_potential=true)
                       push!(energies, kinetic_energy(sys) + U)
                       return nothing
                   end)

        max_drift = maximum(abs.(energies .- E0))
        @test max_drift ≤ 1e-6
    end

    @testset "kinetic energy" begin
        # v = (3,4,0) => v^2 = 25; KE = 0.5 * m * v^2
        positions = [SVector{Dims,Float64}(0.0,0.0,0.0)]
        velocities = [SVector{Dims,Float64}(3.0,4.0,0.0)]
        forces = [SVector{Dims,Float64}(0.0,0.0,0.0)]
        masses = [2.0]
        box = CubicBox(10.0)
        types = [1]
        type_names = Dict(1 => :A)
        sys = System(positions, velocities, forces, masses, box, types, type_names)
        @test isapprox(kinetic_energy(sys), 0.5 * 2.0 * 25.0; atol=1e-12)
    end

    @testset "potential energy conventions" begin
        # A force that supports keyword return_potential
        function f_kw(r; return_potential::Bool=false)
            F = map(x -> -x, r)
            if return_potential

                U = 0.5 * sum(([sum(abs2.(ri)) for ri in r]))  # harmonic potential: 0.5 * Σ r^2
                return (F, U)
            else
                return F
            end
        end
        positions = [SVector{Dims,Float64}(1.0,0.0,0.0)]
        velocities = [SVector{Dims,Float64}(0.0,0.0,0.0)]
        forces = [SVector{Dims,Float64}(0.0,0.0,0.0)]
        masses = [1.0]
        box = CubicBox(10.0)
        types = [1]
        type_names = Dict(1 => :A)
        sys = System(positions, velocities, forces, masses, box, types, type_names)
        @test isapprox(Verlet.Core.potential_energy(sys, f_kw), 0.5; atol=1e-12)

        # A force that does NOT provide potential should error
        f_plain(r) = map(x -> -x, r)
        @test_throws ErrorException Verlet.Core.potential_energy(sys, f_plain)
    end

    @testset "CubicBox + minimum_image!" begin
        box = CubicBox(10.0)
        Δ = [ 6.0, -6.0,  0.1 ]  # in 3D
        Δ = minimum_image(Δ, box)
        @test Δ[1] ==  6.0 - 10.0    # -> -4.0
        @test Δ[2] == -6.0 + 10.0    # ->  4.0
        @test isapprox(Δ[3], 0.1; atol=1e-12)
    end
end

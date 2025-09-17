using Test
using Verlet
using StaticArrays
using LinearAlgebra, StaticArrays

const Dims = 3

struct HarmonicWell end

function Verlet.Core.compute_forces!(::HarmonicWell, sys::System)
    @inbounds for i in 1:natoms(sys)
        sys.forces[i] += -sys.positions[i]
    end
    return sys
end

harmonic_energy(sys::System) = 0.5 * sum(dot(r, r) for r in sys.positions)

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
        positions = [SVector{Dims,Float64}(0.0,0.0,0.0)]
        velocities = [SVector{Dims,Float64}(1.0,0.0,0.0)]
        forces = [SVector{Dims,Float64}(0.0,0.0,0.0)]
        masses = [1.0]
        box = CubicBox(10.0)
        types = [1]
        type_names = Dict(1 => :A)
        sys = System(copy(positions), copy(velocities), copy(forces), copy(masses), box, types, type_names;
                     forcefield=ForceField(()) )
        vv = VelocityVerlet(0.1)
        integrate!(vv, sys, 1)
        @test isapprox(sys.positions[1][1], 0.1; atol=1e-12)
        @test isapprox(sys.velocities[1][1], 1.0; atol=1e-12)

        sys_ho = System([SVector{Dims,Float64}(1.0,0.0,0.0)], [SVector{Dims,Float64}(0.0,0.0,0.0)],
                        [SVector{Dims,Float64}(0.0,0.0,0.0)], [1.0], box, types, type_names;
                        forcefield=ForceField((HarmonicWell(),)))
        integrate!(VelocityVerlet(0.1), sys_ho, 1)
        @test sys_ho.positions[1][1] < 1.0
    end

    @testset "VelocityVerlet integrator behaviour" begin
        function make_system(; pos=SVector{Dims,Float64}(1.0, 0.0, 0.0),
                               vel=SVector{Dims,Float64}(0.0, 1.0, 0.0))
            positions = [pos]
            velocities = [vel]
            forces = [SVector{Dims,Float64}(0.0, 0.0, 0.0)]
            masses = [1.0]
            box = CubicBox(10.0)
            types = [1]
            type_names = Dict(1 => :A)
            return System(copy(positions), copy(velocities), copy(forces), copy(masses), box, types, type_names;
                          forcefield=ForceField((HarmonicWell(),)))
        end
        dt = 0.1

        @testset "matches direct loop" begin
            sys_integrated = make_system()
            sys_reference = deepcopy(sys_integrated)
            vv = VelocityVerlet(dt)
            integrate!(vv, sys_integrated, 5)
            manual_vv = VelocityVerlet(dt)
            for _ in 1:5
                integrate!(manual_vv, sys_reference, 1)
            end
            @test all(isapprox.(sys_integrated.positions, sys_reference.positions; atol=1e-12))
            @test all(isapprox.(sys_integrated.velocities, sys_reference.velocities; atol=1e-12))
        end

        @testset "zero steps no-op" begin
            sys_zero = make_system()
            sys_copy = deepcopy(sys_zero)
            vv = VelocityVerlet(dt)
            calls = Ref(0)
            integrate!(vv, sys_zero, 0; callback = (sys, step) -> (calls[] += 1))
            @test sys_zero.positions == sys_copy.positions
            @test sys_zero.velocities == sys_copy.velocities
            @test calls[] == 0
        end

        @testset "callback early stop" begin
            sys_cb = make_system()
            sys_reference = deepcopy(sys_cb)
            vv = VelocityVerlet(dt)
            last_step = Ref(0)
            integrate!(vv, sys_cb, 10;
                       callback = (sys, step) -> begin
                           last_step[] = step
                           step == 3 ? false : nothing
                       end)
            vv_single = VelocityVerlet(dt)
            for _ in 1:3
                integrate!(vv_single, sys_reference, 1)
            end
            @test last_step[] == 3
            @test all(isapprox.(sys_cb.positions, sys_reference.positions; atol=1e-12))
            @test all(isapprox.(sys_cb.velocities, sys_reference.velocities; atol=1e-12))
        end

        @testset "rejects negative steps" begin
            sys_neg = make_system()
            vv = VelocityVerlet(dt)
            @test_throws ArgumentError integrate!(vv, sys_neg, -2)
        end
    end

    @testset "Supports arbitrary spatial dimension" begin
        D2 = 2
        positions2 = [SVector{D2,Float64}(0.25, -0.4), SVector{D2,Float64}(-1.1, 0.9)]
        velocities2 = [SVector{D2,Float64}(0.3, 0.1), SVector{D2,Float64}(-0.2, 0.05)]
        forces2 = [SVector{D2,Float64}(0.0, 0.0) for _ in 1:2]
        masses2 = [1.0, 2.0]
        box2 = CubicBox(4.0)
        types2 = [1, 1]
        type_names2 = Dict(1 => :A)
        sys2 = System(positions2, velocities2, forces2, masses2, box2, types2, type_names2;
                      forcefield=ForceField((HarmonicWell(),)))
        wrap_positions!(sys2.positions, sys2.box)
        integrate!(VelocityVerlet(0.02), sys2, 10)
        @test length(sys2.positions[1]) == D2
        @test all(isfinite, first(sys2.positions))

        D4 = 4
        positions4 = [SVector{D4,Float64}(ntuple(i -> 0.1 * i, D4)...)]
        velocities4 = [zero(SVector{D4,Float64})]
        forces4 = [zero(SVector{D4,Float64})]
        masses4 = [1.0]
        box4 = CubicBox(6.0)
        types4 = [1]
        type_names4 = Dict(1 => :A)
        sys4 = System(positions4, velocities4, forces4, masses4, box4, types4, type_names4;
                      forcefield=ForceField((HarmonicWell(),)))
        wrap_positions!(sys4.positions, sys4.box)
        integrate!(VelocityVerlet(0.01), sys4, 1)
        @test length(sys4.positions[1]) == D4
    end

    @testset "Conjugate-gradient minimisation" begin
        positions = [SVector{Dims,Float64}(2.0, -1.0, 0.5), SVector{Dims,Float64}(-1.5, 0.7, 1.2)]
        init_positions = copy(positions)
        velocities = [SVector{Dims,Float64}(0.0, 0.0, 0.0) for _ in positions]
        forces = [SVector{Dims,Float64}(0.0, 0.0, 0.0) for _ in positions]
        masses = ones(Float64, length(positions))
        box = CubicBox(20.0)
        types = ones(Int, length(positions))
        type_names = Dict(1 => :A)
        sys = System(copy(positions), copy(velocities), copy(forces), copy(masses), box, types, type_names;
                     forcefield=ForceField((HarmonicWell(),)))

        energies = Float64[]
        energy0 = harmonic_energy(sys)
        cg = ConjugateGradient(harmonic_energy; tol=1e-10)
        integrate!(cg, sys, 200; callback = (s, iter, E) -> push!(energies, E))
        energy_final = harmonic_energy(sys)

        @test energy_final ≤ energy0 - 1e-8
        @test maximum(norm.(sys.positions)) ≤ 1e-5
        @test all(diff(energies) .≤ 1e-10)

        sys2 = System(copy(init_positions), copy(velocities), copy(forces), copy(masses), box, types, type_names)
        @test_throws ArgumentError integrate!(cg, sys2, 10)
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
        sys = System(positions, velocities, forces, masses, box, types, type_names;
                     forcefield=ForceField((HarmonicWell(),)))

        E0 = kinetic_energy(sys) + harmonic_energy(sys)

        energies = Float64[E0]
        dt = 0.002
        steps = 2000
        integrate!(VelocityVerlet(dt), sys, steps;
                   callback = (sys, step) -> begin
                       push!(energies, kinetic_energy(sys) + harmonic_energy(sys))
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

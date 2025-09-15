using Test
using Verlet
using LinearAlgebra
using StaticArrays

function setup_test_system()
    N = 10
    D = 3
    L = 20.0
    box = CubicBox(L)
    positions = [SVector{D,T_Float}(rand(D) .* L .- (L/2)...) for _ in 1:N]
    velocities = [SVector{D,T_Float}(randn(D)...) for _ in 1:N]
    forces = [zero(SVector{D,T_Float}) for _ in 1:N]
    masses = ones(T_Float, N)
    types = ones(T_Int, N)
    type_names = Dict(1 => :A)
    sys = System(positions, velocities, forces, masses, box, types, type_names)
    return sys
end

@testset "LennardJones ForceField" begin
    sys = setup_test_system()

    ϵ = 1.0
    σ = 1.0
    rc = 2.5
    lj_pair = LJPair(ϵ, σ, rc)
    params = PairTable(fill(lj_pair, (1, 1)))
    exclusions = Tuple{T_Int,T_Int}[]
    lj = LennardJones(params, exclusions, 0.5)

    ff = Verlet.Neighbors.ForceField((lj,))

    master_skin = 0.5
    master_nl = MasterNeighborList(master_skin)

    for method in [:cells, :bruteforce, :all_pairs]
        Verlet.Neighbors.build_all_neighbors!(master_nl, ff, sys, method=method)

        @test master_nl isa MasterNeighborList
        @test lj.neighbors isa PotentialNeighborList

        Verlet.Core.compute_all_forces!(sys, ff)

        F_ref = zeros(SVector{3,T_Float}, sys.natoms)
        for i in 1:sys.natoms-1
            for j in i+1:sys.natoms
                Δ = sys.positions[i] - sys.positions[j]
                Δ = minimum_image(Δ, sys.box)
                r2 = dot(Δ, Δ)
                if r2 < rc^2
                    invr2 = 1 / r2
                    s2 = σ^2 * invr2
                    s6 = s2^3
                    fr_over_r = 24 * ϵ * (2 * s6^2 - s6) * invr2
                    fvec = fr_over_r .* Δ
                    F_ref[i] += fvec
                    F_ref[j] -= fvec
                end
            end
        end

        @test sys.forces ≈ F_ref atol=1e-10
    end
end

function setup_test_system_with_charges()
    N = 10
    D = 3
    L = 20.0
    box = CubicBox(L)
    positions = [SVector{D,T_Float}(rand(D) .* L .- (L/2)...) for _ in 1:N]
    velocities = [SVector{D,T_Float}(randn(D)...) for _ in 1:N]
    forces = [zero(SVector{D,T_Float}) for _ in 1:N]
    masses = ones(T_Float, N)
    types = ones(T_Int, N)
    type_names = Dict(1 => :A)
    sys = System(positions, velocities, forces, masses, box, types, type_names)
    return sys
end

@testset "Coulomb ForceField" begin
    sys = setup_test_system_with_charges()

    q1q2 = 1.0
    rc = 5.0
    coul_pair = CoulPair(q1q2, rc)
    params = PairTable(fill(coul_pair, (1, 1)))
    exclusions = Tuple{T_Int,T_Int}[]
    coul = Coulomb(params, exclusions, 0.5)

    ff = Verlet.Neighbors.ForceField((coul,))

    master_skin = 0.5
    master_nl = MasterNeighborList(master_skin)

    for method in [:cells, :bruteforce, :all_pairs]
        Verlet.Neighbors.build_all_neighbors!(master_nl, ff, sys, method=method)
        Verlet.Core.compute_all_forces!(sys, ff)

        F_ref = zeros(SVector{3,T_Float}, sys.natoms)
        for i in 1:sys.natoms-1
            for j in i+1:sys.natoms
                Δ = sys.positions[i] - sys.positions[j]
                Δ = minimum_image(Δ, sys.box)
                r2 = dot(Δ, Δ)
                if r2 < rc^2 && r2 > 0
                    r = sqrt(r2)
                    fr = q1q2 / r2
                    fvec = (fr / r) .* Δ
                    F_ref[i] += fvec
                    F_ref[j] -= fvec
                end
            end
        end

        @test sys.forces ≈ F_ref atol=1e-10
    end
end

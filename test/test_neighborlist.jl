using Test, LinearAlgebra, Verlet, Random, StaticArrays
Random.seed!(0xC0FFEE)

function setup_test_system(N, D, L)
    positions = [SVector{D,Float64}(rand(D) .* L .- (L/2)...) for _ in 1:N]
    velocities = [zero(SVector{D,Float64}) for _ in 1:N]
    forces = [zero(SVector{D,Float64}) for _ in 1:N]
    masses = ones(Float64, N)
    box = CubicBox(L)
    types = ones(Int, N)
    type_names = Dict(1 => :A)
    sys = System(positions, velocities, forces, masses, box, types, type_names)
    wrap_positions!(sys.positions, sys.box)
    return sys
end

@testset "NeighborList" begin
    @testset "MasterNeighborList construction and invariants" begin
        N, D = 10, 3
        L = 8.0
        sys = setup_test_system(N, D, L)

        r_verlet = 2.5 + 0.4
        master_nl = MasterNeighborList(0.4)
        Verlet.Neighbors.build_master_neighborlist!(master_nl, sys.positions, sys.box; r_verlet=r_verlet, method=:cells)

        # All neighbor pairs within r_verlet
        for entry in master_nl.entries
            i, j, r2 = entry.i, entry.j, entry.r2
            Δ = sys.positions[i] - sys.positions[j]
            Δ = minimum_image(Δ, sys.box)
            @test dot(Δ,Δ) ≈ r2
        end
    end

    @testset "LJ forces with NeighborList match O(N²) reference" begin
        N, D = 24, 3
        L = 12.0
        sys = setup_test_system(N, D, L)

        ϵ, σ, rcut, skin = 1.0, 1.0, 2.5, 0.4
        lj_pair = LJPair(ϵ, σ, rcut)
        params = PairTable(fill(lj_pair, (1, 1)))
        exclusions = Tuple{Int,Int}[]
        lj = LennardJones(params, exclusions, skin)
        ff = Verlet.Neighbors.ForceField((lj,))
        master_nl = MasterNeighborList(skin)

        Verlet.Neighbors.build_all_neighbors!(master_nl, ff, sys, method=:cells)
        Verlet.Core.compute_all_forces!(sys, ff)
        F_nl = sys.forces

        U_nl = 0.0
        for pair_info in lj.neighbors
            r2 = Verlet.Neighbors._squared_distance_min_image(sys.positions, pair_info.i, pair_info.j, sys.box)
            if r2 < rcut^2
                inv_r2 = 1 / r2
                s2 = σ^2 * inv_r2
                s6 = s2^3
                U_nl += 4ϵ * (s6^2 - s6)
            end
        end


        # Reference (brute-force)
        F_ref = zeros(SVector{3,Float64}, sys.natoms)
        U_ref = 0.0
        for i in 1:sys.natoms-1
            for j in i+1:sys.natoms
                Δ = sys.positions[i] - sys.positions[j]
                Δ = minimum_image(Δ, sys.box)
                r2 = dot(Δ, Δ)
                if r2 < rcut^2
                    inv_r2 = 1 / r2
                    s2 = σ^2 * inv_r2
                    s6 = s2^3
                    U_ref += 4ϵ * (s6^2 - s6)
                    f_r = 24ϵ * (2s6^2 - s6) * inv_r2
                    f_vec = f_r .* Δ
                    F_ref[i] += f_vec
                    F_ref[j] -= f_vec
                end
            end
        end

        @test isapprox(F_nl, F_ref; atol=1e-10, rtol=1e-10)
        @test isapprox(U_nl,  U_ref; atol=1e-10, rtol=1e-10)
    end
end
using Test, LinearAlgebra, Verlet, Random
Random.seed!(0xC0FFEE)
@testset "NeighborList" begin
    @testset "NeighborList construction and invariants" begin
        N, D = 10, 3
        L = 8.0
        box = CubicBox(L)
        R = randn(N, D)
        wrap_positions!(R, box)

        cutoff = 2.5
        skin   = 0.4
        nlist = build_neighborlist(R, box; cutoff=cutoff, skin=skin)

        # Offsets form valid CSR
        @test length(nlist.offsets) == N + 1
        @test first(nlist.offsets) == 1
        @test last(nlist.offsets) == length(nlist.pairs) + 1

        # All neighbor pairs within rlist
        rlist2 = (cutoff + skin)^2
        for i in 1:N
            for idx in nlist.offsets[i]:(nlist.offsets[i+1]-1)
                j = nlist.pairs[idx]
                Δ = R[i,:] .- R[j,:]
                minimum_image!(Δ, box)
                @test dot(Δ,Δ) <= rlist2 + 1e-12
            end
        end
    end

    @testset "maybe_rebuild! triggers on half-skin" begin
        N, D = 5, 3
        L = 10.0
        box = CubicBox(L)
        R = zeros(N, D)
        cutoff, skin = 2.0, 0.6
        nlist = build_neighborlist(R, box; cutoff=cutoff, skin=skin)

        # Move one particle by < skin/2 : no rebuild
        R[1,1] += 0.25  # < 0.3
        @test maybe_rebuild!(nlist, R, box) == false

        # Cross skin/2 : triggers rebuild
        R[1,1] += 0.1  # now 0.35 > 0.3
        @test maybe_rebuild!(nlist, R, box) == true
    end

    @testset "LJ forces with NeighborList match O(N²) reference" begin
        N, D = 24, 3
        L = 12.0
        box = CubicBox(L)
        R = rand(N, D) .* (L/2) .- (L/4)
        wrap_positions!(R, box)

        ϵ, σ, cutoff, skin = 1.0, 1.0, 2.5, 0.4
        # Build list with rlist = cutoff + skin
        nlist = build_neighborlist(R, box; cutoff=cutoff, skin=skin)

        # Reference (brute-force) and list-aware comparisons
        F_ref, U_ref = lj_forces(R, box; ϵ=ϵ, σ=σ, rcut=cutoff, return_potential=true)
        F_nl,  U_nl  = lj_forces(R, box, nlist; ϵ=ϵ, σ=σ, shift=false, return_potential=true)

        @test isapprox(F_nl, F_ref; atol=1e-10, rtol=1e-10)
        @test isapprox(U_nl,  U_ref; atol=1e-10, rtol=1e-10)
    end

    @testset "NeighborList stays valid under small motion; rebuild changes content" begin
        N, D = 16, 3
        L = 10.0
        box = CubicBox(L)
        R = rand(N, D) .* L .- (L/2)
        cutoff, skin = 2.2, 0.6
        nlist = build_neighborlist(R, box; cutoff=cutoff, skin=skin)

        pairs_before = copy(nlist.pairs)
        offsets_before = copy(nlist.offsets)

        # Small random motion below half-skin → no rebuild
        R .+= 0.1 .* randn(N, D)
        wrap_positions!(R, box)
        @test maybe_rebuild!(nlist, R, box) == false
        @test nlist.pairs == pairs_before
        @test nlist.offsets == offsets_before

        # Larger motion → rebuild likely alters neighbor graph
        R .+= 0.4 .* randn(N, D) # push beyond skin/2
        wrap_positions!(R, box)
        did = maybe_rebuild!(nlist, R, box)
        @test did == true
    end

    @testset "wrap_positions! maps coords to (-L/2, L/2]" begin
        box = CubicBox(5.0)
        R = [ 3.1 -2.6; -4.9 4.9 ]
        wrap_positions!(R, box)
        @test all(-box.L/2 .< R) && all(R .<= box.L/2)
    end

    @testset "Energy shift at cutoff is preserved with NeighborList" begin
        L = 30.0; box = CubicBox(L)
        R = [0.0 0.0 0.0; 2.0 0.0 0.0]
        cutoff, skin = 2.5, 0.4
        nlist = build_neighborlist(R, box; cutoff=cutoff, skin=skin)
        F1, U1 = lj_forces(R, box, nlist; rcut=cutoff, shift=true, return_potential=true)
        R2 = [0.0 0.0 0.0; 2.51 0.0 0.0]
        nlist2 = build_neighborlist(R2, box; cutoff=cutoff, skin=skin)
        F2, U2 = lj_forces(R2, box, nlist2; rcut=cutoff, shift=true, return_potential=true)
        @test isapprox(norm(F2), 0.0; atol=1e-12)
        @test U2 ≈ 0.0 atol=1e-10
        @test U1 < 0.0
    end
end
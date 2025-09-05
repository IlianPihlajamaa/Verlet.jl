using Test
using LinearAlgebra
using Verlet
using Random

@testset "CellGrid" begin
    @testset "CellGrid build and rebin invariants" begin
        rng = MersenneTwister(0xC311C311)
        N, D = 100, 3
        L = 12.0
        box = CubicBox(L)
        R = rand(rng, N, D) .* L .- (L/2)
        wrap_positions!(R, box)

        rcut, skin = 2.5, 0.4
        rlist = rcut + skin

        grid = build_cellgrid(R, box; cell_size=rlist)
        @test grid.L == L
        @test all(x -> x ≥ 1, grid.dims)
        @test length(grid.next) == N
        @test length(grid.heads) == prod(grid.dims)

        rebin!(grid, R, box)
        @test length(grid.next) == N
    end

    @testset "Cell-based neighbor list ≡ naive rlist filter (half list)" begin
        rng = MersenneTwister(0xB15B15B1)
        N, D = 128, 3
        L = 15.0
        box = CubicBox(L)
        R = rand(rng, N, D) .* L .- (L/2)
        wrap_positions!(R, box)

        cutoff, skin = 2.2, 0.5
        nlist = build_neighborlist_cells(R, box; cutoff=cutoff, skin=skin)

        rlist2 = (cutoff + skin)^2
        ref_pairs = Vector{Tuple{Int,Int}}()
        for i in 1:N-1, j in i+1:N
            # Create a proper mutable vector for minimum_image! without @view on a broadcast
            Δ = Vector{Float64}(undef, 3)
            @inbounds begin
                Δ[1] = R[i,1] - R[j,1]
                Δ[2] = R[i,2] - R[j,2]
                Δ[3] = R[i,3] - R[j,3]
            end
            minimum_image!(Δ, box)
            if dot(Δ,Δ) ≤ rlist2 + 1e-12
                push!(ref_pairs, (i,j))
            end
        end
        sort!(ref_pairs)

        got = Vector{Tuple{Int,Int}}()
        for i in 1:N
            for idx in nlist.offsets[i]:(nlist.offsets[i+1]-1)
                j = nlist.pairs[idx]
                @test j > i
                push!(got, (i,j))
            end
        end
        sort!(got)
        @test got == ref_pairs
    end

    @testset "LJ with half list from cells matches brute force" begin
        rng = MersenneTwister(0xF4AF1157)
        N, D = 64, 3
        L = 12.0
        box = CubicBox(L)
        R = rand(rng, N, D) .* L .- (L/2)
        wrap_positions!(R, box)

        ϵ, σ, cutoff, skin = 1.0, 1.0, 2.5, 0.4
        nlist = build_neighborlist_cells(R, box; cutoff=cutoff, skin=skin)

        F_ref, U_ref = lj_forces(R, box; ϵ=ϵ, σ=σ, rcut=cutoff, return_potential=true)
        F_nl,  U_nl  = lj_forces(R, box, nlist; ϵ=ϵ, σ=σ, shift=false, rcut=cutoff, return_potential=true)

        @test isapprox(F_nl, F_ref; atol=1e-10, rtol=1e-10)
        @test isapprox(U_nl,  U_ref; atol=1e-10, rtol=1e-10)
    end

    @testset "Rebin + maybe_rebuild! integration" begin
        rng = MersenneTwister(0xC3C1C3B1)
        N, D = 80, 3
        L = 14.0
        box = CubicBox(L)
        R = rand(rng, N, D) .* L .- (L/2)
        wrap_positions!(R, box)

        cutoff, skin = 2.5, 0.4
        grid = build_cellgrid(R, box; cell_size=cutoff+skin)
        nlist = build_neighborlist_cells(R, box; cutoff=cutoff, skin=skin, grid=grid)

    # Small bounded move: ensure max per-particle displacement < skin/2
    δ = 0.9 * (skin/2) / sqrt(3)   # 10% safety margin under half-skin
    R .+= δ .* (2 .* rand(rng, N, D) .- 1)
    wrap_positions!(R, box)
    @test maybe_rebuild!(nlist, R, box) == false

        R .+= 0.35 .* randn(rng, N, D); wrap_positions!(R, box)
        if maybe_rebuild!(nlist, R, box)
            rebin!(grid, R, box)
            nlist2 = build_neighborlist_cells(R, box; cutoff=cutoff, skin=skin, grid=grid)
            @test length(nlist2.pairs) > 0
        end
    end
end
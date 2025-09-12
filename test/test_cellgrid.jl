using Test
using LinearAlgebra, StaticArrays
using Verlet
using Random

@testset "CellGrid" begin
    @testset "CellGrid build and rebin invariants" begin
        rng = MersenneTwister(0xC311C311)
        N, D = 100, 3
        L = 12.0
        box = CubicBox(L)
        R = rand(rng, SVector{D, Float64}, N) .* L 
        for i in 1:N
            R[i] -= (L/2) * ones(SVector{D, Float64})
        end
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

    @testset "Cell-based master neighbor list" begin
        rng = MersenneTwister(0xB15B15B1)
        N, D = 128, 3
        L = 15.0
        box = CubicBox(L)
        R = rand(rng, SVector{D, Float64}, N) .* L 
        for i in 1:N
            R[i] -= (L/2) * ones(SVector{D, Float64})
        end
        wrap_positions!(R, box)

        cutoff, skin = 2.2, 0.5
        r_verlet = cutoff
        master_nl = MasterNeighborList(skin)
        build_master_neighborlist!(master_nl, R, box; r_verlet=r_verlet, method=:cells)

        rlist2 = (cutoff + skin)^2
        ref_pairs = Tuple{Int,Int}[]
        for i in 1:N-1, j in i+1:N
            Δ = R[i] - R[j]
            Δ = minimum_image(Δ, box)
            if dot(Δ,Δ) ≤ rlist2
                push!(ref_pairs, (i,j))
            end
        end
        sort!(ref_pairs)

        got_pairs = [(p.i, p.j) for p in master_nl.entries]
        sort!(got_pairs)

        @test got_pairs == ref_pairs
    end
end
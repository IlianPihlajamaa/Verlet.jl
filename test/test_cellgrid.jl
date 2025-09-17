using Test
using LinearAlgebra, StaticArrays
using Verlet
using Random

@testset "CellGrid" begin
    @testset "CellGrid build and rebin invariants" begin
        for D in (1, 2, 3, 4)
            @testset "D=$D" begin
                rng = MersenneTwister(0xC311C311)
                N = 100
                L = 12.0
                box = CubicBox(L)
                R = rand(rng, SVector{D, Float64}, N) .* L
                half = (L / 2) * ones(SVector{D, Float64})
                @inbounds for i in 1:N
                    R[i] -= half
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
        end
    end

    @testset "Cell-based master neighbor list" begin
        for D in (1, 2, 3, 4)
            @testset "D=$D" begin
                rng = MersenneTwister(0xB15B15B1)
                N = 128
                L = 15.0
                box = CubicBox(L)
                R = rand(rng, SVector{D, Float64}, N) .* L
                half = (L / 2) * ones(SVector{D, Float64})
                @inbounds for i in 1:N
                    R[i] -= half
                end
                wrap_positions!(R, box)

                cutoff, skin = 2.2, 0.5
                r_verlet = cutoff + skin
                master_nl = MasterNeighborList(R, box; cutoff=cutoff, skin=skin)
                build_master_neighborlist!(master_nl, R, box; r_verlet=r_verlet, method=:cells)

                rlist2 = (cutoff + skin)^2
                ref_pairs = Tuple{Int,Int}[]
                for i in 1:N-1, j in i+1:N
                    Δ = R[i] - R[j]
                    Δ = minimum_image(Δ, box)
                    if dot(Δ, Δ) ≤ rlist2
                        push!(ref_pairs, (i, j))
                    end
                end
                sort!(ref_pairs)

                got_pairs = [Tuple(p) for p in master_nl.pairs]
                sort!(got_pairs)

                @test got_pairs == ref_pairs
            end
        end
    end
end

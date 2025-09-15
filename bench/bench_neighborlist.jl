# Benchmark neighbor-list construction vs. force evaluation
#
# Usage:
#   julia bench/bench_neighborlist.jl
#   julia bench/bench_neighborlist.jl -- N=128,256,512 D=3 rho=1.0 rcut=2.5 skin=0.4 samples=10 builder=both potential=false
#
# The script activates the local bench environment (BenchmarkTools).

import Pkg; Pkg.activate(@__DIR__)

using BenchmarkTools
using Random
using LinearAlgebra, StaticArrays
using Verlet
using Printf

# ----------------------------
# CLI options (very lightweight)
# ----------------------------
function parse_args(args::Vector{String})
    opts = Dict{String,Any}(
        "N" => [128, 256, 512, 1024, 2048, 4096, 8192],
        "D" => 3,
        "rho" => 1.0,
        "rcut" => 1.25,
        "skin" => 0.4,
        "samples" => 10,
        "seed" => 42,
        "potential" => false,      # set true to include potential accumulation
        "builder" => "all",       # "all_pairs" | "cells" | "brute_force" | "all"
        "storage" => "both",   # "entries" | "csr" | "both"
    )
    for arg in args
        if occursin("=", arg)
            k, v = split(arg, "=", limit=2)
            if k == "N"
                opts["N"] = parse.(Int, split(v, ","))
            elseif k == "D"
                opts["D"] = parse(Int, v)
            elseif k == "rho"
                opts["rho"] = parse(Float64, v)
            elseif k == "rcut"
                opts["rcut"] = parse(Float64, v)
            elseif k == "skin"
                opts["skin"] = parse(Float64, v)
            elseif k == "samples"
                opts["samples"] = parse(Int, v)
            elseif k == "seed"
                opts["seed"] = parse(Int, v)
            elseif k == "potential"
                opts["potential"] = lowercase(v) in ("1","true","yes","y")
            elseif k == "builder"
                opts["builder"] = lowercase(v)
            elseif k == "storage"
                opts["storage"] = lowercase(v)
            end
        end
    end
    return opts
end

const OPTS = parse_args(copy(ARGS))
Random.seed!(OPTS["seed"])

# ----------------------------
# Helpers
# ----------------------------
make_positions(N::Int, D::Int, L::Float64) = L*rand(SVector{D, Float64}, N) .- (L/2*ones(SVector{D, Float64}),)

ns(t) = minimum(t).time

function fmt_ns(x)
    if x < 1_000
        return "$(x) ns"
    elseif x < 1_000_000
        return @sprintf("%.3f μs", x/1_000)
    elseif x < 1_000_000_000
        return @sprintf("%.3f ms", x/1_000_000)
    else
        return @sprintf("%.3f s", x/1_000_000_000)
    end
end

function setup_system(N, D, L)
    positions = make_positions(N, D, L)
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

# ----------------------------
# Bench cases
# ----------------------------
function bench_case(sys, ff, master_nl; method::Symbol, samples::Int)
    build_trial = @benchmark Verlet.Neighbors.build_all_neighbors!($master_nl, $ff, $sys, method=$method) samples=samples evals=1
    Verlet.Neighbors.build_all_neighbors!(master_nl, ff, sys, method=method)

    force_trial = @benchmark Verlet.Core.compute_all_forces!($sys, $ff) samples=samples evals=1

    lj = ff.layers[1]
    avg_deg = length(lj.neighborlist.neighbors) / sys.natoms

    return (; build_trial, force_trial, avg_deg)
end

# ----------------------------
# Run
# ----------------------------
println("Benchmarking Neighbor Lists")
println("Params: D=$(OPTS["D"]), rho=$(OPTS["rho"]), rcut=$(OPTS["rcut"]), skin=$(OPTS["skin"]), samples=$(OPTS["samples"]), potential=$(OPTS["potential"]), builder=$(OPTS["builder"]), storage=$(OPTS["storage"])\n")

for N in OPTS["N"]
    D   = OPTS["D"]
    rho = OPTS["rho"]
    L   = (N / rho)^(1/D)

    sys = setup_system(N, D, L)

    # Setup ForceField
    ϵ, σ, rcut, skin = 1.0, 1.0, OPTS["rcut"], OPTS["skin"]
    lj_pair = LJPair(ϵ, σ, rcut)
    params = PairTable(fill(lj_pair, (1, 1)))
    exclusions = Tuple{Int,Int}[]
    lj = LennardJones(params, exclusions, skin)
    ff = Verlet.Neighbors.ForceField((lj,))
    storages = OPTS["storage"] == "both" ? ("entries", "csr") : (OPTS["storage"],)
    nl =  Verlet.Neighbors.MasterNeighborCSRList(skin; N=sys.natoms)
    Verlet.Neighbors.build_all_neighbors!(nl, ff, sys, method=:cells)
    @profview for i = 1:100
        Verlet.Neighbors.build_all_neighbors!(nl, ff, sys, method=:cells)
        Verlet.Core.compute_all_forces!(sys, ff)
    end


    for storage in storages
        master_nl = storage == "csr" ? Verlet.Neighbors.MasterNeighborCSRList(skin; N=sys.natoms) : Verlet.Neighbors.MasterNeighborList(skin)

        if OPTS["builder"] in ("all_pairs", "all")
            res = bench_case(sys, ff, master_nl; method=:all_pairs, samples=OPTS["samples"])
            t_build = ns(res.build_trial)
            t_force = ns(res.force_trial)

            println("=== N=$N (all_pairs, $storage) ===")
            println(rpad("avg_deg", 16), rpad("build", 14), rpad("force(NL)", 14))
            println(rpad(@sprintf("%.2f", res.avg_deg), 16),
                    rpad(fmt_ns(t_build), 14),
                    rpad(fmt_ns(t_force), 14))
            println()
        end

        if OPTS["builder"] in ("cells", "all")
            res = bench_case(sys, ff, master_nl; method=:cells, samples=OPTS["samples"])
            t_build = ns(res.build_trial)
            t_force = ns(res.force_trial)

            println("=== N=$N (cells, $storage) ===")
            println(rpad("avg_deg", 16), rpad("build", 14), rpad("force(NL)", 14))
            println(rpad(@sprintf("%.2f", res.avg_deg), 16),
                    rpad(fmt_ns(t_build), 14),
                    rpad(fmt_ns(t_force), 14))
            println()
        end

        if OPTS["builder"] in ("brute_force", "all")
            res = bench_case(sys, ff, master_nl; method=:bruteforce, samples=OPTS["samples"])
            t_build = ns(res.build_trial)
            t_force = ns(res.force_trial)

            println("=== N=$N (bruteforce, $storage) ===")
            println(rpad("avg_deg", 16), rpad("build", 14), rpad("force(NL)", 14))
            println(rpad(@sprintf("%.2f", res.avg_deg), 16),
                    rpad(fmt_ns(t_build), 14),
                    rpad(fmt_ns(t_force), 14))
            println()
        end
    end
end

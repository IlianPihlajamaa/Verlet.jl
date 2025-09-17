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
        "N" => 10 .^ (2:6),  # list of system sizes
        "D" => 3,
        "rho" => 1.0,
        "rcut" => 1.25,
        "skin" => 0.4,
        "samples" => 10,
        "seed" => 42,
        "potential" => false,      # set true to include potential accumulation
        "builder" => "cells",       # "all_pairs" | "cells" | "brute_force" | "all"
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
println("Params: D=$(OPTS["D"]), rho=$(OPTS["rho"]), rcut=$(OPTS["rcut"]), skin=$(OPTS["skin"]), samples=$(OPTS["samples"]), potential=$(OPTS["potential"]), builder=$(OPTS["builder"])\n")




println(rpad("N", 22), rpad("avg_deg", 16), rpad("build", 14), rpad("force(NL)", 14))

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
    master_nl = Verlet.Neighbors.MasterNeighborList(sys; cutoff=rcut, skin=skin)

    Verlet.Neighbors.build_all_neighbors!(master_nl, ff, sys; method=:cells)
    Verlet.Core.compute_all_forces!(sys, ff)  # warmup

    if OPTS["builder"] in ("all_pairs", "all")
        res = bench_case(sys, ff, master_nl; method=:all_pairs, samples=OPTS["samples"])
        t_build = ns(res.build_trial)
        t_force = ns(res.force_trial)

        println(rpad("N=$N (all_pairs)", 22),
                rpad(@sprintf("%.2f", res.avg_deg), 16),
                rpad(fmt_ns(t_build), 14),
                rpad(fmt_ns(t_force), 14))
    end

    if OPTS["builder"] in ("cells", "all")
        res = bench_case(sys, ff, master_nl; method=:cells, samples=OPTS["samples"])
        t_build = ns(res.build_trial)
        t_force = ns(res.force_trial)

        println(rpad("N=$N (cells)", 22),
                rpad(@sprintf("%.2f", res.avg_deg), 16),
                rpad(fmt_ns(t_build), 14),
                rpad(fmt_ns(t_force), 14))
    end

    if OPTS["builder"] in ("brute_force", "all")
        res = bench_case(sys, ff, master_nl; method=:bruteforce, samples=OPTS["samples"])
        t_build = ns(res.build_trial)
        t_force = ns(res.force_trial)

        println(rpad("N=$N (bruteforce)", 22),
                rpad(@sprintf("%.2f", res.avg_deg), 16),
                rpad(fmt_ns(t_build), 14),
                rpad(fmt_ns(t_force), 14))
    end
end


## laptop
# Params: D=3, rho=1.0, rcut=1.25, skin=0.4, samples=10, potential=false, builder=cells

# N                     avg_deg         build         force(NL)     
# N=100 (cells)         9.33            59.600 μs     9.100 μs      
# N=1000 (cells)        9.44            764.000 μs    144.900 μs    
# N=10000 (cells)       9.38            7.847 ms      1.698 ms      
# N=100000 (cells)      9.40            90.758 ms     19.574 ms     
# N=1000000 (cells)     9.41            1.495 s       372.285 ms

## phys-computee006

# Params: D=3, rho=1.0, rcut=1.25, skin=0.4, samples=10, potential=false, builder=cells

# N                     avg_deg         build         force(NL)     
# N=100 (cells)         9.33            94.245 μs     13.264 μs     
# N=1000 (cells)        9.44            1.326 ms      199.467 μs    
# N=10000 (cells)       9.38            13.805 ms     2.136 ms      
# N=100000 (cells)      9.40            158.474 ms    22.870 ms     
# N=1000000 (cells)     9.41            1.776 s       264.127 ms
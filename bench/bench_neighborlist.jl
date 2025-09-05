# Benchmark neighbor-list construction vs. force evaluation
#
# Usage:
#   julia bench/bench_neighborlist.jl
#   julia bench/bench_neighborlist.jl -- N=128,256,512 D=3 L=20 rcut=2.5 skin=0.4 samples=10
#
# The script activates the local bench environment (BenchmarkTools).


using BenchmarkTools
using Random
using LinearAlgebra
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
        "rcut" => 2.5,
        "skin" => 0.4,
        "samples" => 10,
        "seed" => 42,
        "potential" => false, # set true to include potential accumulation
    )
    for arg in args
        if occursin("=", arg)
            k, v = split(arg, "=", limit=2)
            if k == "N"
                opts["N"] = parse.(Int, split(v, ","))
            elseif k == "D"
                opts["D"] = parse(Int, v)
            elseif k == "L"
                opts["L"] = parse(Float64, v)
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
            end
        end
    end
    return opts
end

const OPTS = parse_args(copy(ARGS))
Random.seed!(OPTS["seed"])

function make_positions(N::Int, D::Int, L::Float64)
    # Uniform in (-L/2, L/2]
    R = (rand(N, D) .* L) .- (L/2)
    return R
end

function bench_case(N::Int; D::Int, rho::Float64, rcut::Float64, skin::Float64, samples::Int, wantU::Bool)
    L = (N / rho)^(1/D)
    box = CubicBox(L)
    R = make_positions(N, D, L)
    wrap_positions!(R, box)

    # --- Build neighbor list (measured) ---
    build_trial = @benchmark build_neighborlist($R, $box; cutoff=$rcut, skin=$skin) samples=samples evals=1
    nlist = build_neighborlist(R, box; cutoff=rcut, skin=skin)  # realize once for force benchmarks

    # neighbor stats (average degree)
    n_pairs = length(nlist.pairs)
    avg_deg = n_pairs / N

    # --- Force eval with neighbor list (measured) ---
    if wantU
        force_nl_trial = @benchmark lj_forces($R, $box, $nlist; rcut=$rcut, return_potential=true) samples=samples evals=1
    else
        force_nl_trial = @benchmark lj_forces($R, $box, $nlist; rcut=$rcut, return_potential=false) samples=samples evals=1
    end

    # --- Brute-force reference (measured) ---
    if wantU
        force_bf_trial = @benchmark lj_forces($R, $box; rcut=$rcut, return_potential=true) samples=samples evals=1
    else
        force_bf_trial = @benchmark lj_forces($R, $box; rcut=$rcut, return_potential=false) samples=samples evals=1
    end

    return (; build_trial, force_nl_trial, force_bf_trial, avg_deg)
end

function ns(t)
    return minimum(t).time
end

function fmt_ns(x::Integer)
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

println("Benchmarking Neighbor List (construction vs. force eval)")
println("Params: D=$(OPTS["D"]), rho=$(OPTS["rho"]), rcut=$(OPTS["rcut"]), skin=$(OPTS["skin"]), samples=$(OPTS["samples"]), potential=$(OPTS["potential"])")
println()
println(rpad("N", 8), rpad("avg_deg", 12), rpad("build", 14), rpad("force(NL)", 14), rpad("force(O(N^2))", 16), "speedup (bf/nl)")
println("-"^80)

for N in OPTS["N"]
    res = bench_case(N; D=OPTS["D"], rho=OPTS["rho"], rcut=OPTS["rcut"], skin=OPTS["skin"],
                     samples=OPTS["samples"], wantU=OPTS["potential"])
    t_build  = ns(res.build_trial)
    t_nl     = ns(res.force_nl_trial)
    t_brutef = ns(res.force_bf_trial)
    speedup  = t_brutef / t_nl
    println(rpad(string(N), 8),
            rpad(@sprintf("%.2f", res.avg_deg), 12),
            rpad(fmt_ns(t_build), 14),
            rpad(fmt_ns(t_nl), 14),
            rpad(fmt_ns(t_brutef), 16),
            @sprintf("%.2f×", speedup))
end

println()
println("Notes:")
println("• 'avg_deg' ≈ average neighbor count per particle at rcut+skin (symmetric list).")
println("• 'force(NL)' measures lj_forces(R, box, nlist; ...), excluding build time.")
println("• 'force(O(N^2))' is the brute-force reference over all pairs.")
println("• Times are minimum over 'samples' runs (BenchmarkTools).")
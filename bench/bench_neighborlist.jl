# Benchmark neighbor-list construction vs. force evaluation
#
# Usage:
#   julia bench/bench_neighborlist.jl
#   julia bench/bench_neighborlist.jl -- N=128,256,512 D=3 rho=1.0 rcut=2.5 skin=0.4 samples=10 builder=both potential=false
#
# The script activates the local bench environment (BenchmarkTools).

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
        "rcut" => 2.5,
        "skin" => 0.4,
        "samples" => 10,
        "seed" => 42,
        "potential" => false,      # set true to include potential accumulation
        "builder" => "both",       # "naive" | "cells" | "both"
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
make_positions(N::Int, D::Int, L::Float64) = (rand(N, D) .* L) .- (L/2)

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

# ----------------------------
# Bench cases
# ----------------------------
function bench_case_naive(R, box; rcut, skin, samples, wantU::Bool)
    build_trial = @benchmark build_neighborlist($R, $box; cutoff=$rcut, skin=$skin) samples=samples evals=1
    nlist = build_neighborlist(R, box; cutoff=rcut, skin=skin)

    # average degree for symmetric CSR list
    avg_deg = length(nlist.pairs) / size(R,1)

    if wantU
        force_trial = @benchmark lj_forces($R, $box, $nlist; rcut=$rcut, return_potential=true, shift=false) samples=samples evals=1
    else
        force_trial = @benchmark lj_forces($R, $box, $nlist; rcut=$rcut, return_potential=false, shift=false) samples=samples evals=1
    end
    return (; build_trial, force_trial, avg_deg)
end

function bench_case_cells(R, box; rcut, skin, samples, wantU::Bool)
    cell_size = rcut + skin

    # Grid build (binning)
    grid_build_trial = @benchmark build_cellgrid($R, $box; cell_size=$cell_size) samples=samples evals=1
    grid = build_cellgrid(R, box; cell_size=cell_size)

    # Half-list build using an existing grid
    nl_from_grid_trial = @benchmark build_neighborlist_cells($R, $box; cutoff=$rcut, skin=$skin, grid=$grid) samples=samples evals=1
    nlc = build_neighborlist_cells(R, box; cutoff=rcut, skin=skin, grid=grid)

    # Half-list build letting the builder manage the grid internally (total cost)
    nl_nogrid_trial = @benchmark build_neighborlist_cells($R, $box; cutoff=$rcut, skin=$skin) samples=samples evals=1

    # average degree (report symmetric equivalent for comparability)
    avg_deg_sym = (2 * length(nlc.pairs)) / size(R,1)

    if wantU
        force_trial = @benchmark lj_forces($R, $box, $nlc; rcut=$rcut, return_potential=true, shift=false) samples=samples evals=1
    else
        force_trial = @benchmark lj_forces($R, $box, $nlc; rcut=$rcut, return_potential=false, shift=false) samples=samples evals=1
    end

    return (; grid_build_trial, nl_from_grid_trial, nl_nogrid_trial, force_trial, avg_deg_sym)
end

function bench_case_bf(R, box; rcut, samples, wantU::Bool)
    if wantU
        return @benchmark lj_forces($R, $box; rcut=$rcut, return_potential=true, shift=false) samples=samples evals=1
    else
        return @benchmark lj_forces($R, $box; rcut=$rcut, return_potential=false, shift=false) samples=samples evals=1
    end
end

# ----------------------------
# Run
# ----------------------------
println("Benchmarking Neighbor Lists (CSR vs. cell-linked half-list)")
println("Params: D=$(OPTS["D"]), rho=$(OPTS["rho"]), rcut=$(OPTS["rcut"]), skin=$(OPTS["skin"]), samples=$(OPTS["samples"]), potential=$(OPTS["potential"]), builder=$(OPTS["builder"])\n")

for N in OPTS["N"]
    D   = OPTS["D"]
    rho = OPTS["rho"]
    L   = (N / rho)^(1/D)
    box = CubicBox(L)
    R   = make_positions(N, D, L)
    wrap_positions!(R, box)

    # Always get brute-force for speedup comparison
    bf_trial = bench_case_bf(R, box; rcut=OPTS["rcut"], samples=OPTS["samples"], wantU=OPTS["potential"])
    t_bf = ns(bf_trial)

    if OPTS["builder"] in ("naive", "both")
        resN = bench_case_naive(R, box; rcut=OPTS["rcut"], skin=OPTS["skin"], samples=OPTS["samples"], wantU=OPTS["potential"])
        t_build = ns(resN.build_trial)
        t_force = ns(resN.force_trial)
        speedup = t_bf / t_force

        println("=== N=$N (CSR naive) ===")
        println(rpad("avg_deg", 16), rpad("build", 14), rpad("force(NL)", 14), rpad("force(O(N^2))", 16), "speedup (bf/nl)")
        println(rpad(@sprintf("%.2f", resN.avg_deg), 16),
                rpad(fmt_ns(t_build), 14),
                rpad(fmt_ns(t_force), 14),
                rpad(fmt_ns(t_bf), 16),
                @sprintf("%.2f×", speedup))
        println()
    end

    if OPTS["builder"] in ("cells", "both")
        resC = bench_case_cells(R, box; rcut=OPTS["rcut"], skin=OPTS["skin"], samples=OPTS["samples"], wantU=OPTS["potential"])
        t_grid   = ns(resC.grid_build_trial)
        t_nl_g   = ns(resC.nl_from_grid_trial)
        t_nl_all = ns(resC.nl_nogrid_trial)
        t_force  = ns(resC.force_trial)
        speedup  = t_bf / t_force

        println("=== N=$N (Cell-linked half-list) ===")
        println(rpad("avg_deg(sym)", 16), rpad("build(cells+nl)", 18), rpad("grid", 14), rpad("nl_from_grid", 16),
                rpad("force(NL)", 14), rpad("force(O(N^2))", 16), "speedup (bf/nl)")
        println(rpad(@sprintf("%.2f", resC.avg_deg_sym), 16),
                rpad(fmt_ns(t_nl_all), 18),
                rpad(fmt_ns(t_grid), 14),
                rpad(fmt_ns(t_nl_g), 16),
                rpad(fmt_ns(t_force), 14),
                rpad(fmt_ns(t_bf), 16),
                @sprintf("%.2f×", speedup))
        println()
    end
end

println("Notes:")
println("• CSR 'avg_deg' is the true symmetric average degree at (rcut+skin).")
println("• Cell 'avg_deg(sym)' is 2×(half-list pairs)/N to make it comparable to CSR.")
println("• 'build(cells+nl)' runs build_neighborlist_cells without supplying a grid (includes grid build).")
println("• 'grid' and 'nl_from_grid' show the components when reusing a prebuilt CellGrid.")
println("• 'force(NL)' always uses shift=false here to match the brute-force settings.")


# === N=128 (CSR naive) ===
# avg_deg         build         force(NL)     force(O(N^2))   speedup (bf/nl)
# 94.62           1.098 ms      299.387 μs    398.961 μs      1.33×

# === N=128 (Cell-linked half-list) ===
# avg_deg(sym)    build(cells+nl)   grid          nl_from_grid    force(NL)     force(O(N^2))   speedup (bf/nl)
# 94.62           236.830 μs        3.875 μs      221.405 μs      121.154 μs    398.961 μs      3.29×

# === N=256 (CSR naive) ===
# avg_deg         build         force(NL)     force(O(N^2))   speedup (bf/nl)
# 101.30          3.050 ms      609.003 μs    1.515 ms        2.49×

# === N=256 (Cell-linked half-list) ===
# avg_deg(sym)    build(cells+nl)   grid          nl_from_grid    force(NL)     force(O(N^2))   speedup (bf/nl)
# 241.52          2.447 ms          7.983 μs      2.511 ms        626.953 μs    1.515 ms        2.42×

# === N=512 (CSR naive) ===
# avg_deg         build         force(NL)     force(O(N^2))   speedup (bf/nl)
# 102.77          8.699 ms      1.212 ms      5.568 ms        4.59×

# === N=512 (Cell-linked half-list) ===
# avg_deg(sym)    build(cells+nl)   grid          nl_from_grid    force(NL)     force(O(N^2))   speedup (bf/nl)
# 206.79          7.532 ms          15.869 μs     7.641 ms        1.161 ms      5.568 ms        4.80×

# === N=1024 (CSR naive) ===
# avg_deg         build         force(NL)     force(O(N^2))   speedup (bf/nl)
# 102.27          28.956 ms     2.287 ms      21.251 ms       9.29×

# === N=1024 (Cell-linked half-list) ===
# avg_deg(sym)    build(cells+nl)   grid          nl_from_grid    force(NL)     force(O(N^2))   speedup (bf/nl)
# 102.27          8.866 ms          30.799 μs     8.958 ms        1.271 ms      21.251 ms       16.73×

# === N=2048 (CSR naive) ===
# avg_deg         build         force(NL)     force(O(N^2))   speedup (bf/nl)
# 101.92          103.745 ms    4.563 ms      85.157 ms       18.66×

# === N=2048 (Cell-linked half-list) ===
# avg_deg(sym)    build(cells+nl)   grid          nl_from_grid    force(NL)     force(O(N^2))   speedup (bf/nl)
# 101.92          17.241 ms         61.166 μs     17.444 ms       2.590 ms      85.157 ms       32.88×

# === N=4096 (CSR naive) ===
# avg_deg         build         force(NL)     force(O(N^2))   speedup (bf/nl)
# 102.35          373.021 ms    9.059 ms      339.137 ms      37.44×

# === N=4096 (Cell-linked half-list) ===
# avg_deg(sym)    build(cells+nl)   grid          nl_from_grid    force(NL)     force(O(N^2))   speedup (bf/nl)
# 102.35          36.573 ms         122.675 μs    37.585 ms       5.318 ms      339.137 ms      63.77×

# === N=8192 (CSR naive) ===
# avg_deg         build         force(NL)     force(O(N^2))   speedup (bf/nl)
# 102.07          1.416 s       17.979 ms     1.348 s         74.96×

# === N=8192 (Cell-linked half-list) ===
# avg_deg(sym)    build(cells+nl)   grid          nl_from_grid    force(NL)     force(O(N^2))   speedup (bf/nl)
# 102.07          86.460 ms         245.081 μs    87.178 ms       10.648 ms     1.348 s         126.56×
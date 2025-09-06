# Verlet Neighbor List & Potential API 

## 1. Goals

* **Modular, composable potentials** (LJ, Coulomb, etc.).
* Each potential owns a **type-stable, contiguous neighbor list** .
* Single **master neighbor list** for efficient spatial search.
* Exclusions and per-pair cutoffs applied during potential neighbor list construction.
* r² is computed on-the-fly; optional stale r² in `masterNL` used only for coarse filtering. Valid only at construction.
* Master NL specifies skin, and computes **maximum total Verlet radius**.
* Memory-efficient, vectorizable, and fully type-stable kernels.

---

## 2. Core Data Structures

### 2.1 Master Neighbor List

* Built using the **largest Verlet radius** across all potentials:

$$
r_\text{Verlet}^\text{master} = \max_\text{pot} \max_\text{pairs}(rc_\text{pair}) + \text{skin}_\text{master}
$$

* Stores candidate pairs `(i,j)` and optionally stale `r²` for early filtering:

```julia
struct MasterNeighborEntry
    i::Int
    j::Int
    r2::Float64
end

struct MasterNeighborList{FloatT<:Number}
    skin::FloatT
    entries::Vector{MasterNeighborEntry}
    etc
end
```

* Built **once per rebuild**, typically using spatial partitioning or cell lists. CSR storage.

---

### 2.2 Per-Potential Neighbor List

* Owned by the potential; type-stable and contiguous StructArrays:

```julia
struct NeighborPair{F<:AbstractPotentialPair, IntT<:Integer}
    i::IntT
    j::IntT
    pair::F  # concrete pair parameters including rc
end
masterNL = Vector{NeighborPair{<:AbstractPotentialPair}}()
```

* Each `pair` struct includes **all force parameters and the per-pair cutoff**:

```julia
struct LJPair{FloatT<:Number}
    ε::FloatT
    σ::FloatT
    rc::FloatT
end
```

* Filtering from `masterNL` to `potNL` uses something like:

```julia
function build_neighbors_from_master!(pot::Potential, sys::System, state::State, masterNL)
    ivec, jvec, pairvec = Int[], Int[], F[] # Use preallocated buffers at rebuild time. empty! -> push!
    for entry in masterNL
        p = pot.params[sys.types[entry.i], sys.types[entry.j]]
        if !is_excluded(pot, entry.i, entry.j) && entry.r2 < (p.rc + pot.skin)^2
            push!(ivec, entry.i)
            push!(jvec, entry.j)
            push!(pairvec, p)
        end
    end
    pot.neighbors = PotentialNeighborList(ivec, jvec, pairvec) # StructArrays
end
```

* **Notes:**

  * `r²` in masterNL is used only for coarse filtering.
  * potNL is **contiguous, type-stable**, and ready for vectorized kernels.

---

### 2.3  Pair Potential & Pair API

```julia
abstract type AbstractPairPotential end

struct LennardJones{IntT<:Integer, T<:AbstractPotentialPair} <: AbstractPairPotential
    params::PairTable{T}        # per-pair parameters including rc
    exclusions::Vector{Tuple{IntT, IntT}}
    neighbors::PotentialNeighborList{T}
end

struct Coulomb <: AbstractPairPotential
    params::PairTable{CoulPair}
    skin::Float64
    neighbors::PotentialNeighborList{CoulPair}
end
```

* `PairTable` is typically a dictionary keyed by `(type_i, type_j)`:

```julia
struct PairTable{F<:AbstractPotentialPair}
    table::Matrix{F}
end
```

---

### 2.4 ForceField & Simulation API

```julia
struct ForceField{ForcesTuple}
    layers::ForcesTuple  # type-stable tuple containing different forces (lj, coul, etc)
end

function build_all_neighbors!(ff::ForceField, sys::System)
    # 1. Compute master NL radius
    rc_max = maximum(maximum([p.rc for p in values(pot.params.table)]) + pot.skin
                      for pot in ff.layers)
    masterNL = build_master_neighborlist(sys, rc_max)
    
    # 2. Build per-potential NLs
    for pot in ff.layers
        build_neighbors_from_master!(pot, sys, masterNL)
    end
    
    return masterNL
end

function compute_all_forces!(ff::ForceField, state::State)
    for pot in ff.layers
        compute_forces!(pot, state)  # kernel operates on pot.neighbors
    end
end
```

---

## 3. Neighbor List Rebuild Logic

* **Master NL** rebuild triggered when **any particle moves > skin / 2**.
* **Per-potential NL** filtered from master:

```julia
for pot in ff.layers
    build_neighbors_from_master!(pot, sys, masterNL)
end
```

* Efficient: only one spatial search (master NL), per-potential lists are cheap linear filtering.

---

### 5. Usage Example

```julia
sys = System(...)  # atoms, positions, types

lj = LennardJones(params=lj_pairs)
coul = Coulomb(params=coul_pairs)

ff = ForceField(layers=(lj, coul))

masterNL = build_all_neighbors!(ff, sys)
compute_all_forces!(ff, state)
```

* Master NL uses **max(rc + skin)**.
* Per-potential NLs are filtered and type-stable.
* Kernels operate on contiguous arrays; r² computed on-the-fly.

---

✅ **Outcome**

* Modular, type-stable neighbor lists per potential.
* Single master NL reduces rebuild cost.
* Per-pair cutoffs and skins fully supported.
* StructArrays enable efficient CPU/GPU execution.
* Exclusions handled early; inner loops are branchless and vectorizable.

#


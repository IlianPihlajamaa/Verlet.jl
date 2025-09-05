Awesome — here’s a complete documentation pass that (a) upgrades your docstrings and (b) gives you a minimal-but-solid Documenter.jl site with runnable examples, performance notes, and pitfalls.

---

# ✅ Patch `src/Verlet.jl` (docstrings only)

Paste the following over your current file (only docstrings/typos changed; code logic untouched). This also adds a module-level docstring.

````julia
"""
    module Verlet

Basic velocity Verlet integrator for tiny molecular-dynamics style problems.

# Overview

- **Data layout:** Positions and velocities are `N×D` `Matrix{Float64}` (rows = particles, columns = spatial dimensions).
- **Masses:** Length-`N` `Vector{Float64}`.
- **Forces:** User-supplied function that takes positions `(N×D)` and returns forces `(N×D)` in the same units.

# Units

You pick units, but they must be **consistent**:
- Positions `r` in length (e.g., m, Å),
- Velocities `v` in length/time,
- Masses `m` in mass,
- Forces `F` in mass·length/time²,
- Timestep `dt` in time.

# Stability & Pitfalls

- Choose `dt` small enough relative to the fastest vibrational frequency in your system.
- Use `Float64` throughout to minimize drift.
- Energy is approximately conserved for small `dt`, but will drift if `dt` is too large.

# Exports

[`ParticleSystem`](@ref), [`velocity_verlet!`](@ref), [`kinetic_energy`](@ref), [`potential_energy`](@ref)
"""
module Verlet

export ParticleSystem, velocity_verlet!, kinetic_energy, potential_energy

"""
    ParticleSystem

A lightweight container for particle states.

# Fields
- `positions::Matrix{Float64}`: `(N × D)` array of positions (rows = particles).
- `velocities::Matrix{Float64}`: `(N × D)` array of velocities.
- `masses::Vector{Float64}`: length-`N` vector of particle masses.

!!! tip "Memory layout"
    For best performance, keep `positions` and `velocities` as `Matrix{Float64}` with fixed `D` (1–3 typical).
"""
mutable struct ParticleSystem
    positions::Matrix{Float64}   # (N × D)
    velocities::Matrix{Float64}  # (N × D)
    masses::Vector{Float64}      # (N)
end

"""
    velocity_verlet!(system::ParticleSystem, forces::Function, dt::Float64)

Advance `system` by one timestep `dt` using the velocity Verlet integrator.

The `forces` function **must** accept a positions matrix `(N × D)` and return an
`(N × D)` matrix of forces (same shape).

!!! info "Algorithm"
    Given positions `r`, velocities `v`, accelerations `a = F(r)/m`, and `dt`:
    1. `r ← r + v*dt + 0.5*a*dt^2`
    2. Recompute `a' = F(r)/m`
    3. `v ← v + 0.5*(a + a')*dt`

# Examples
```julia
julia> using Verlet

julia> forces(r) = zeros(size(r))  # free particle
forces (generic function with 1 method)

julia> ps = ParticleSystem([0.0 0.0], [1.0 0.0], [1.0]);

julia> velocity_verlet!(ps, forces, 0.1);  # advances in-place

julia> ps.positions[1,1]
0.1
````

!!! warning "Shape rules"
\- `size(positions) == size(velocities)`
\- `length(masses) == size(positions, 1)`
"""
function velocity\_verlet!(system::ParticleSystem, forces::Function, dt::Float64)
@assert size(system.positions) == size(system.velocities) "positions/velocities must be same size"
@assert length(system.masses) == size(system.positions, 1) "length(masses) must equal number of particles"

```
dt2 = dt * dt

# Current accelerations a = F/m (broadcast per-row / per-particle)
a = forces(system.positions) ./ system.masses

# Update positions
system.positions .+= system.velocities .* dt .+ 0.5 .* a .* dt2

# New accelerations from updated positions
a_new = forces(system.positions) ./ system.masses

# Update velocities
system.velocities .+= 0.5 .* (a .+ a_new) .* dt

return system
```

end

"""
kinetic\_energy(system::ParticleSystem) -> Float64

Total kinetic energy: `∑ ½ mᵢ ‖vᵢ‖²`.

# Examples

```julia
julia> ps = ParticleSystem([0.0 0.0; 1.0 0.0],
                           [1.0 2.0; 0.0 0.0],
                           [1.0, 2.0]);

julia> kinetic_energy(ps)
0.5*(1*(1^2+2^2) + 2*(0^2+0^2))  # = 2.5
2.5
```

"""
function kinetic\_energy(system::ParticleSystem)::Float64
@assert size(system.positions) == size(system.velocities) "positions/velocities must be same size"
@assert length(system.masses) == size(system.positions, 1) "length(masses) must equal number of particles"
v2 = sum(abs2, system.velocities; dims=2)
return 0.5 \* sum(system.masses .\* vec(v2))
end

"""
potential\_energy(system::ParticleSystem, forces::Function) -> Float64

Try to obtain **total potential energy** from the user-supplied `forces` function.

Supported convention:

* If `forces(r; return_potential=true)` is supported, it must return `(F, U)`
  where `F` is the `(N × D)` force matrix and `U` is a scalar total potential energy.

If this convention is not supported, an error is thrown.

# Examples

```julia
julia> function ho_forces(r; return_potential=false)
           F = -r
           return return_potential ? (F, 0.5 * sum(abs2, r)) : F
       end;

julia> ps = ParticleSystem([1.0 0.0], [0.0 0.0], [1.0]);

julia> potential_energy(ps, ho_forces)
0.5
```

"""
function potential\_energy(system::ParticleSystem, forces::Function)::Float64
F\_U\_kw = try
forces(system.positions; return\_potential=true)
catch
nothing
end
if F\_U\_kw !== nothing && F\_U\_kw isa Tuple && length(F\_U\_kw) == 2
\_, U = F\_U\_kw
return float(U)
end
error("Force function does not support `return_potential=true`; cannot compute potential\_energy.")
end

end # module Verlet

````

---

# 📖 New documentation site (Documenter.jl)

Create these files:

## `docs/Project.toml`
```toml
[deps]
Documenter = "e30172f5-a6a5-5a46-863b-614d45cd2de4"

[compat]
Documenter = "1"
````

## `docs/make.jl`

```julia
using Documenter
using Verlet

DocMeta.setdocmeta!(Verlet, :DocTestSetup, :(using Verlet); recursive=true)

makedocs(
    sitename = "Verlet.jl",
    modules = [Verlet],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", "false") == "true"),
    pages = [
        "Home" => "index.md",
        "Guide" => Any[
            "Forces & Potentials" => "guide/forces.md",
            "Numerics & Pitfalls" => "guide/numerics.md",
        ],
        "API" => "api.md",
    ],
)

deploydocs(; repo = "github.com/yourname/Verlet.jl", devbranch = "main")
```

## `docs/src/index.md`

````markdown
# Verlet.jl

A minimal **velocity Verlet** integrator for tiny MD-style problems.

!!! note
    You choose the units; just keep them consistent across positions `r`, velocities `v`, masses `m`,
    forces `F`, and timestep `dt`.

## Quickstart

```@example quickstart
using Verlet

# Free particle in 2D
forces(r) = zeros(size(r))

ps = ParticleSystem([0.0 0.0],  # 1×2 positions
                    [1.0 0.0],  # 1×2 velocities
                    [1.0])      # masses

dt = 0.1
velocity_verlet!(ps, forces, dt)
ps.positions
````

## Harmonic oscillator

```@example ho
using Verlet

# Hooke's law with k = 1, potential U = 0.5 * |r|^2
function ho_forces(r; return_potential=false)
    F = -r
    return return_potential ? (F, 0.5 * sum(abs2, r)) : F
end

ps = ParticleSystem([1.0 0.0], [0.0 0.0], [1.0])

dt = 0.1
for _ in 1:100
    velocity_verlet!(ps, ho_forces, dt)
end

(kin = kinetic_energy(ps), pot = potential_energy(ps, ho_forces))
```

## Energy monitoring

```@example energy
using Verlet

function ho_forces(r; return_potential=false)
    F = -r
    return return_potential ? (F, 0.5 * sum(abs2, r)) : F
end

ps = ParticleSystem([1.0 0.0], [0.0 1.0], [1.0])
dt = 0.05

energies = Float64[]
for _ in 1:200
    velocity_verlet!(ps, ho_forces, dt)
    push!(energies, kinetic_energy(ps) + potential_energy(ps, ho_forces))
end

(round(minimum(energies), digits=6), round(maximum(energies), digits=6))
```

## Performance tips

* Keep arrays as `Matrix{Float64}` / `Vector{Float64}` to avoid type instability.
* Prefer **in-place** force computations in your own code paths; if you must allocate, reuse buffers.
* Avoid huge `dt`. Start small (e.g., `1e-3` in your time units) and increase cautiously.

See also: \[Numerics & Pitfalls]\(@ref numerics).

````

## `docs/src/api.md`
```markdown
# API Reference

```@docs
Verlet
ParticleSystem
velocity_verlet!
kinetic_energy
potential_energy
````

````

## `docs/src/guide/forces.md`
```markdown
# Forces & Potentials

A **force function** takes positions `r::Matrix{Float64}` of shape `(N×D)` and returns a force matrix of the same shape.

```@example forces
# Free particle (no forces)
forces_free(r) = zeros(size(r))

# Linear spring to the origin (Hooke's law, k = 1)
forces_ho(r) = -r

# With potential-energy support:
function forces_ho_with_U(r; return_potential=false)
    F = -r
    return return_potential ? (F, 0.5 * sum(abs2, r)) : F
end
nothing
````

To make [`potential_energy`](@ref) work, support the keyword protocol:

```@example forces
# Called internally as: forces(r; return_potential=true) => (F, U)
nothing
```

If your force does **not** provide `return_potential=true`, `potential_energy` will throw with a clear message.

````

## `docs/src/guide/numerics.md`
```markdown
# Numerics & Pitfalls  {#numerics}

## Timestep stability

- Velocity Verlet is symplectic and stable for sufficiently small `dt`.
- If you see exploding energies or trajectories, reduce `dt`.

## Precision

Use `Float64` across positions, velocities, and masses to limit drift.

## Shapes & broadcasting

- Arrays are `(N×D)` with rows as particles.
- `masses` is length `N`. Broadcasting `forces(r) ./ masses` divides each row by its particle mass.

## Boundary conditions (not included)

No periodic or reflective boundaries are included. Implement them in your force or position update logic if needed.

## Potential energy protocol

- Provide `forces(r; return_potential=true) => (F, U)` to enable `potential_energy`.
- If unknown, prefer returning only `F` and skip potential-energy reporting.
````

---

# 🧪 Docs build & tests (optional)

Once you’ve created those files, you can run:

* Run tests:

```bash
julia -e 'import Pkg; Pkg.test()'
```

* Build docs locally:

```bash
julia --project=docs -e 'import Pkg; Pkg.activate("docs"); Pkg.instantiate(); include("docs/make.jl")'
```

---

If you want, I can also add a tiny `README.md` snippet showing badges and the Quickstart section mirrored from the docs.

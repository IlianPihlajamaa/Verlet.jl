"""

* langevin_baoab!(ps::ParticleSystem, forces, dt; γ, T, kB=1.0, rng=Random.default_rng())
*

Advance the system by **one NVT step** using the **Langevin BAOAB** integrator.
 
This scheme integrates the Langevin SDE
 
`
m dv = [F(r) - γ m v] dt + √(2γ m kB T) dW,     dr = v dt
`
 
via the **B → A → O → A → B** splitting with an **exact Ornstein–Uhlenbeck (O)** step.
It is widely used for robust thermostatting and near-optimal configurational sampling.
When `γ == 0`, the O-step is the identity and BAOAB reduces to velocity–Verlet (up to roundoff).
 
# Arguments
- `ps::ParticleSystem`: System with `positions::Vector{SVector{Dims,T}}`, `velocities::Vector{SVector{Dims,T}}`, and `masses::Vector{T}`.
- `forces`: A callable `F(R)::Vector{SVector{Dims,T}}` returning forces evaluated at positions `R`.
- `dt::Real`: Time step (same time units as used for `γ`).
 
# Keywords
- `γ::Real`: Friction coefficient **[1/time]**.
- `T::Real`: Target temperature. Must be consistent with `kB` (so that `kB*T` has units of energy).
- `kB::Real = 1.0`: Boltzmann constant in your unit system.
- `rng::AbstractRNG = Random.default_rng()`: RNG used only in the stochastic O-step for reproducibility.
 
# Details
The O-step uses the exact OU update:
`
c = exp(-γ*dt)
sᵢ = sqrt((1 - c^2) * kB*T / mᵢ)  # per-particle, per-component scale
vᵢ[k] ← c * vᵢ[k] + sᵢ * ξ,    where ξ ~ 𝒩(0,1)
`
For numerical stability when `γ*dt ≪ 1`, implementations should use the linearized forms
`c ≈ 1 - γ*dt` and `1 - c^2 ≈ 2γ*dt` to avoid catastrophic cancellation.
 
# Performance tips
- Avoid allocations in the inner loop. Reuse force buffers if you add multi-step drivers.
- Use `@views` for in-place operations on `ps.positions`/`ps.velocities`.
- Pass an explicit `rng` to ensure reproducibility across runs and Julia sessions.
 
# Pitfalls
- **Units matter:** ensure `γ` is in `1/time`, and `kB*T` is in energy, consistent with `r`, `v`, and `m`.
- Very large `γ*dt` yields overdamped dynamics (OK for sampling, not for kinetics).
- With heterogeneous masses, compute `sᵢ` **per particle**; do not share a single scale across species.
- For tiny systems, instantaneous temperature fluctuates strongly — average over time for diagnostics.
 
# Example
```julia
# Pseudocode showing typical usage (assumes your package exports ParticleSystem and BAOAB):
using Random
N, D = 64, 3
ps = ParticleSystem(zeros(N, D), randn(N, D), ones(N))
dt, γ, T = 0.005, 1.0, 1.5
rng = MersenneTwister(123)
forces(R) = -0.1 .* R  # harmonic "bath"
 
for _ in 1:2_000

* langevin_baoab!(ps, forces, dt; γ=γ, T=T, kB=1.0, rng=rng)
  end
  ```
  """
  langevin_baoab!
*

"""

* instantaneous_temperature(ps; kB=1.0) -> Float64
*

Compute the **instantaneous scalar temperature** via equipartition:
 
`
KE = 0.5 * Σᵢ mᵢ ||vᵢ||²
T  = 2 * KE / (kB * dof)
`
 
where `dof = degrees_of_freedom(ps)`. For now, `dof = N * D` (translational DoF only).
 
# Keywords
- `kB::Real = 1.0`: Boltzmann constant in your unit system.
 
# Returns
- `Float64`: instantaneous temperature.
 
# Notes
- For small `N`, expect large fluctuations. Compare **time averages** to the target `T` when thermostatting.
"""
instantaneous_temperature

"""

* degrees_of_freedom(ps) -> Int
*

Return the **effective translational degrees of freedom** of the system.
Currently this is simply `N * D` for `N` particles in `D` spatial dimensions.
 
This function is a hook for future extensions (e.g., constraints, rigid bodies, or COM removal)
that reduce the effective DoF used in temperature calculations.
"""
degrees_of_freedom

"""

* velocity_rescale!(ps, T; kB=1.0)
*

Deterministically rescale all velocities to match a target temperature `T`.
This is a **one-shot utility** often used for quick pre-equilibration or initializing an NVT run.
It is **not** a thermostat and does not generate the correct kinetic energy distribution by itself.
 
# Keywords
- `kB::Real = 1.0`: Boltzmann constant.
 
# Notes
Let `T₀ = instantaneous_temperature(ps; kB)` and `λ = √(T / max(T₀, eps()))`.
This function applies `vᵢ ← λ * vᵢ` in-place.
"""
velocity_rescale!
 # Thermostats and temperature utilities (BAOAB Langevin + helpers)
 #
 # Public API (exported from src/Verlet.jl):
 #   - degrees_of_freedom(ps) -> Int
 #   - instantaneous_temperature(ps; kB=1.0) -> Float64
 #   - velocity_rescale!(ps, T; kB=1.0)
 #   - langevin_baoab!(ps, forces, dt; γ, T, kB=1.0, rng=Random.default_rng())
 #
 # Notes:
 # * Assumes `ps.positions::AbstractMatrix{<:Real}` (N×D),
 #           `ps.velocities::AbstractMatrix{<:Real}` (N×D),
 #           `ps.masses::AbstractVector{<:Real}` (length N).
 # * `forces` is a callable returning an N×D array of forces given positions.

using Random
using ..Constraints: DistanceConstraints, apply_rattle!, apply_shake!

"""
    degrees_of_freedom(ps; constraints=nothing, remove_com=false) -> Int

Effective translational DoF, reduced by number of constraints and optionally COM removal.
This is the canonical method. A no-keyword fallback is provided for backward compatibility.
"""
function degrees_of_freedom(ps; constraints=nothing, remove_com::Bool=false)::Int
    N = length(ps.positions)
    D = length(ps.positions[1])
    dof = N * D
    if constraints !== nothing
        dof -= length(constraints.r0)
    end
    if remove_com
        dof -= D
    end
    return max(dof, 0)
end

"""
    degrees_of_freedom(ps) -> Int

Backward-compatible method: returns `N*D` (no constraints, no COM removal).
"""


"""
    instantaneous_temperature(ps; kB=1.0) -> Float64

Compute instantaneous temperature via equipartition:
T = 2 * KE / (kB * dof).
"""
function instantaneous_temperature(ps; kB::T_float=one(T_float))::T_float where {T_float}
    v = ps.velocities
    m = ps.masses
    KE = T_float(0.5) * sum(m[i] * sum(abs2.(v[i])) for i in eachindex(v))   # scalar
    dof = degrees_of_freedom(ps)
    return (T_float(2.0) * KE) / (kB * dof)
end

"""
    velocity_rescale!(ps, T; kB=1.0)

Deterministically rescale velocities to match target temperature T.
"""
function velocity_rescale!(ps, T::Real; kB::Real=1.0)
    Tinst = instantaneous_temperature(ps; kB=kB)
    λ = sqrt(T / max(Tinst, eps()))
    ps.velocities .*=  λ
    return ps
end

# Safe evaluation of c = exp(-γ*dt) with a short-series fallback
@inline function _ou_coeff(γ::Real, dt::Real)
    x = γ * dt
    if x < 1e-8
        # 1 - x + x^2/2 - x^3/6: accurate near zero and avoids cancellation
        return 1 - x + 0.5x^2 - (x^3)/6
    else
        return exp(-x)
    end
end

"""
    langevin_baoab!(ps, forces, dt; γ, T, kB=1.0, rng=Random.default_rng())

Advance one step with the Langevin BAOAB integrator:
B(half) → A(half) → O(OU) → A(half) → B(half).
With γ=0, this reduces to velocity-Verlet (deterministic).

Arguments
---------
* `ps`     : particle system (positions N×D, velocities N×D, masses length N)
* `forces` : callable F = forces(positions)::N×D
* `dt`     : time step

Keywords
--------
* `γ`  : friction (1/time)
* `T`  : target temperature
* `kB` : Boltzmann constant (default 1.0)
* `rng`: AbstractRNG for reproducibility
"""
function langevin_baoab!(ps, forces, dt; γ, T, kB::Real=1.0, rng::AbstractRNG=Random.default_rng())
    R = ps.positions
    V = ps.velocities
    m = ps.masses
    N = length(R)
    D = length(R[1])

    # Precompute and cache
    invm = 1.0 ./ m                            # length N
    c = _ou_coeff(γ, dt)
    s2_common = max(0.0, 1 - c^2) * (kB * T)   # scalar >= 0

    # Force at start
    F = forces(R)

    # B: half kick   v += (dt/2) * F/m
    for i in 1:N
        V[i] += (0.5 * dt) * F[i] * invm[i]
    end

    # A: half drift  r += (dt/2) * v
    for i in 1:N
        R[i] += (0.5 * dt) * V[i]
    end

    # O: Ornstein–Uhlenbeck (exact)
    if s2_common == 0.0
        # no noise term; pure damping (c==1 ⇒ identity)
        for i in 1:N
            V[i] = c * V[i]
        end
    else
        for i in 1:N
            ξ = randn(rng, SVector{D,Float64}) # N(0,1) for each component
            s = sqrt(s2_common / m[i])
            V[i] = c * V[i] + ξ * s
        end
    end

    # A: half drift
    for i in 1:N
        R[i] += (0.5 * dt) * V[i]
    end

    # Recompute forces at new positions
    F = forces(R)

    # B: half kick
    for i in 1:N
        V[i] += (0.5 * dt) * F[i] * invm[i]
    end

    return ps
end

# ------------------------------------------------------------
# Constrained Langevin (BAOAB with SHAKE/RATTLE projections)
# ------------------------------------------------------------
"""
    langevin_baoab_constrained!(ps::ParticleSystem, forces, dt, cons::DistanceConstraints;
                                γ, T, kB=1.0, rng=Random.default_rng())

Advance one **constrained NVT** BAOAB step with SHAKE/RATTLE projections.

Splitting and projections:

1. B (half kick)         → `apply_rattle!` (velocities)
2. A (half drift)        → `apply_shake!` (positions)
3. O (OU stochastic)     → `apply_rattle!` (velocities)
4. A (half drift)        → `apply_shake!` (positions)
5. Recompute forces
6. B (half kick)         → `apply_rattle!` (velocities)

Notes:
- Uses `exp(-γ*dt)` for the OU decay. Noise variance is `(1 - c^2) * kB*T / m_i` per component.
- With constraints, prefer this method over the unconstrained `langevin_baoab!`.
"""
function langevin_baoab_constrained!(ps::ParticleSystem, forces, dt, cons::DistanceConstraints;
                                     γ, T, kB::Real=1.0, rng::AbstractRNG=Random.default_rng())
    @assert dt > 0
    R = ps.positions
    V = ps.velocities
    m = ps.masses
    N = length(R)
    D = length(R[1])
    # Per-particle inverse masses (broadcast across velocity components)
    invm = 1.0 ./ m

    # Initial forces
    F = forces(R)

    # --- B: half kick
    for i in 1:N
        V[i] += (dt/2) * F[i] * invm[i]
    end
    apply_rattle!(ps, cons)  # enforce Ċ = 0

    # --- A: half drift
    for i in 1:N
        R[i] += (dt/2) * V[i]
    end
    apply_shake!(ps, cons, dt/2)  # enforce C = 0

    # --- O: OU stochastic velocity step
    c = exp(-γ*dt)
    for i in 1:N
        ξ = randn(SVector{D,Float64})  # N(0,1) for each component
        σ = sqrt((1 - c^2) * (kB * T) * invm[i])
        V[i] = c * V[i] + ξ * σ
    end
    apply_rattle!(ps, cons)

    # --- A: half drift
    for i in 1:N
        R[i] += (dt/2) * V[i]
    end
    apply_shake!(ps, cons, dt/2)

    # Recompute forces at new positions
    F = forces(R)

    # --- B: half kick
    for i in 1:N
        V[i] += (dt/2) * F[i] * invm[i]
    end
    apply_rattle!(ps, cons)

    return ps
end

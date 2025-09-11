using Random

"""
    degrees_of_freedom(sys; constraints=nothing, remove_com=false) -> Int

Effective translational DoF, reduced by number of constraints and optionally COM removal.
This is the canonical method. A no-keyword fallback is provided for backward compatibility.
"""
function degrees_of_freedom(sys::System; constraints=nothing, remove_com::Bool=false)::Int
    N = length(sys.positions)
    D = length(sys.positions[1])
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
    instantaneous_temperature(sys; kB=1.0) -> Float64

Compute instantaneous temperature via equipartition:
T = 2 * KE / (kB * dof).
"""
function instantaneous_temperature(sys::System{T,IT,Dims}; kB::T=one(T))::T where {T,IT,Dims}
    v = sys.velocities
    m = sys.masses
    KE = T(0.5) * sum(m[i] * sum(abs2.(v[i])) for i in eachindex(v))   # scalar
    dof = degrees_of_freedom(sys)
    return (T(2.0) * KE) / (kB * dof)
end

"""
    velocity_rescale!(sys, T; kB=1.0)

Deterministically rescale velocities to match target temperature T.
"""
function velocity_rescale!(sys::System, T::Real; kB::Real=1.0)
    Tinst = instantaneous_temperature(sys; kB=kB)
    λ = sqrt(T / max(Tinst, eps()))
    sys.velocities .*=  λ
    return sys
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
    langevin_baoab!(sys, forces, dt; γ, T, kB=1.0, rng=Random.default_rng())

Advance one step with the Langevin BAOAB integrator:
B(half) → A(half) → O(OU) → A(half) → B(half).
With γ=0, this reduces to velocity-Verlet (deterministic).

Arguments
---------
* `sys`     : particle system
* `forces` : callable F = forces(positions)
* `dt`     : time step

Keywords
--------
* `γ`  : friction (1/time)
* `T`  : target temperature
* `kB` : Boltzmann constant (default 1.0)
* `rng`: AbstractRNG for reproducibility
"""
function langevin_baoab!(sys::System{T,IT,Dims}, forces, dt; γ, temp, kB::Real=1.0, rng::AbstractRNG=Random.default_rng()) where {T,IT,Dims}
    R = sys.positions
    V = sys.velocities
    m = sys.masses
    N = length(R)
    D = length(R[1])

    # Precompute and cache
    invm = 1.0 ./ m                            # length N
    c = _ou_coeff(γ, dt)
    s2_common = max(0.0, 1 - c^2) * (kB * temp)   # scalar >= 0

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
            ξ = randn(rng, SVector{D,T}) # N(0,1) for each component
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

    return sys
end

# ------------------------------------------------------------
# Constrained Langevin (BAOAB with SHAKE/RATTLE projections)
# ------------------------------------------------------------
"""
    langevin_baoab_constrained!(sys::System, forces, dt, cons::DistanceConstraints;
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
function langevin_baoab_constrained!(sys::System{T,IT,Dims}, forces, dt, cons::DistanceConstraints;
                                     γ, temp, kB::Real=1.0, rng::AbstractRNG=Random.default_rng()) where {T,IT,Dims}
    @assert dt > 0
    R = sys.positions
    V = sys.velocities
    m = sys.masses
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
    apply_rattle!(sys, cons)  # enforce Ċ = 0

    # --- A: half drift
    for i in 1:N
        R[i] += (dt/2) * V[i]
    end
    apply_shake!(sys, cons, dt/2)  # enforce C = 0

    # --- O: OU stochastic velocity step
    c = exp(-γ*dt)
    for i in 1:N
        ξ = randn(SVector{D,T})  # N(0,1) for each component
        σ = sqrt((1 - c^2) * (kB * temp) * invm[i])
        V[i] = c * V[i] + ξ * σ
    end
    apply_rattle!(sys, cons)

    # --- A: half drift
    for i in 1:N
        R[i] += (dt/2) * V[i]
    end
    apply_shake!(sys, cons, dt/2)

    # Recompute forces at new positions
    F = forces(R)

    # --- B: half kick
    for i in 1:N
        V[i] += (dt/2) * F[i] * invm[i]
    end
    apply_rattle!(sys, cons)

    return sys
end

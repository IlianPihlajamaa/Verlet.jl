"""
    LangevinBAOAB(dt; γ, temp, kB=1.0, wrap=false, rng=Random.default_rng())

Stochastic BAOAB Langevin integrator (unconstrained). Uses the `System`'s
`ForceField` for force evaluations.
"""
struct LangevinBAOAB{T,RNG} <: AbstractIntegrator
    dt::T
    gamma::T
    temp::T
    kB::T
    wrap::Bool
    rng::RNG
end

function LangevinBAOAB(dt::Real; γ::Real, temp::Real, kB::Real=1.0,
                       wrap::Bool=false, rng::AbstractRNG=Random.default_rng())
    T = promote_type(typeof(dt), typeof(γ), typeof(temp), typeof(kB))
    LangevinBAOAB{T,typeof(rng)}(T(dt), T(γ), T(temp), T(kB), wrap, rng)
end

"""
    LangevinBAOABConstrained(dt, constraints; γ, temp, kB=1.0, wrap=false, rng=Random.default_rng())

Constrained Langevin BAOAB integrator using SHAKE/RATTLE projections.
"""
struct LangevinBAOABConstrained{T,RNG} <: AbstractIntegrator
    dt::T
    gamma::T
    temp::T
    kB::T
    wrap::Bool
    rng::RNG
    constraints::DistanceConstraints
end

function LangevinBAOABConstrained(dt::Real, constraints::DistanceConstraints;
                                  γ::Real, temp::Real, kB::Real=1.0,
                                  wrap::Bool=false, rng::AbstractRNG=Random.default_rng())
    T = promote_type(typeof(dt), typeof(γ), typeof(temp), typeof(kB))
    LangevinBAOABConstrained{T,typeof(rng)}(T(dt), T(γ), T(temp), T(kB), wrap, rng, constraints)
end

@inline function _ou_coeff(γ::Real, dt::Real)
    x = γ * dt
    if x < 1e-8
        return 1 - x + 0.5x^2 - (x^3)/6
    else
        return exp(-x)
    end
end

function step!(integrator::LangevinBAOAB, system::System)
    dt = integrator.dt
    γ = integrator.gamma
    temp = integrator.temp
    kB = integrator.kB
    wrap = integrator.wrap
    rng = integrator.rng

    R = system.positions
    V = system.velocities
    m = system.masses
    N = length(R)
    D = length(R[1])
    invm = 1.0 ./ m

    forces = _update_forces!(system)

    for i in 1:N
        V[i] += (0.5 * dt) * forces[i] * invm[i]
    end

    for i in 1:N
        R[i] += (0.5 * dt) * V[i]
    end
    wrap && wrap_positions!(R, system.box)

    c = _ou_coeff(γ, dt)
    s2_common = max(zero(kB), 1 - c^2) * (kB * temp)
    if s2_common == 0
        for i in 1:N
            V[i] = c * V[i]
        end
    else
        for i in 1:N
            ξ = randn(rng, SVector{D, typeof(R[i][1])})
            s = sqrt(s2_common / m[i])
            V[i] = c * V[i] + ξ * s
        end
    end

    for i in 1:N
        R[i] += (0.5 * dt) * V[i]
    end
    wrap && wrap_positions!(R, system.box)

    forces = _update_forces!(system)
    for i in 1:N
        V[i] += (0.5 * dt) * forces[i] * invm[i]
    end

    return nothing
end

function step!(integrator::LangevinBAOABConstrained, system::System)
    dt = integrator.dt
    γ = integrator.gamma
    temp = integrator.temp
    kB = integrator.kB
    wrap = integrator.wrap
    rng = integrator.rng
    cons = integrator.constraints

    R = system.positions
    V = system.velocities
    m = system.masses
    N = length(R)
    D = length(R[1])
    invm = 1.0 ./ m

    forces = _update_forces!(system)

    for i in 1:N
        V[i] += (0.5 * dt) * forces[i] * invm[i]
    end
    apply_rattle!(system, cons)

    for i in 1:N
        R[i] += (0.5 * dt) * V[i]
    end
    apply_shake!(system, cons, dt/2)
    wrap && wrap_positions!(R, system.box)

    c = exp(-γ * dt)
    for i in 1:N
        ξ = randn(rng, SVector{D, typeof(R[i][1])})
        σ = sqrt(max(zero(kB), (1 - c^2) * kB * temp * invm[i]))
        V[i] = c * V[i] + ξ * σ
    end
    apply_rattle!(system, cons)

    for i in 1:N
        R[i] += (0.5 * dt) * V[i]
    end
    apply_shake!(system, cons, dt/2)
    wrap && wrap_positions!(R, system.box)

    forces = _update_forces!(system)
    for i in 1:N
        V[i] += (0.5 * dt) * forces[i] * invm[i]
    end
    apply_rattle!(system, cons)

    return nothing
end

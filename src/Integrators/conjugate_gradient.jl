mutable struct CGState{V,T}
    grad::Vector{V}
    dir::Vector{V}
    backup::Vector{V}
    grad_tmp::Vector{V}
    energy::T
    initialized::Bool
    stop::Bool
end

"""
    ConjugateGradient(energy; tol=1e-8, alpha0=1.0, min_alpha=1e-8, c1=1e-4, wrap=true)

Polak–Ribière non-linear conjugate-gradient minimiser. Forces are obtained from
the system's `ForceField`; the user supplies an `energy(system)` callable used
for the Armijo backtracking line search. Each call to `step!` performs one CG
iteration; the integrator maintains internal state so subsequent steps continue
from where the previous iteration stopped.
"""
mutable struct ConjugateGradient{T,E} <: AbstractIntegrator
    tol::T
    alpha0::T
    min_alpha::T
    c1::T
    wrap::Bool
    energy::E
    state::CGState
end

function ConjugateGradient(energy::Function; tol::Real=1e-8, alpha0::Real=1.0,
                           min_alpha::Real=1e-8, c1::Real=1e-4, wrap::Bool=true)
    T = promote_type(typeof(tol), typeof(alpha0), typeof(min_alpha), typeof(c1))
    dummy_state = CGState(Vector{Any}(), Vector{Any}(), Vector{Any}(), Vector{Any}(), T(0), false, false)
    ConjugateGradient{T,typeof(energy)}(T(tol), T(alpha0), T(min_alpha), T(c1), wrap, energy, dummy_state)
end

function _ensure_state!(integrator::ConjugateGradient, system::System)
    V = typeof(system.positions[1])
    Tval = typeof(integrator.tol)
    state = integrator.state
    if !(state isa CGState{V,Tval}) || length(state.grad) != length(system.positions)
        integrator.state = CGState{V,Tval}(
            [zero(V) for _ in system.positions],
            [zero(V) for _ in system.positions],
            copy(system.positions),
            [zero(V) for _ in system.positions],
            zero(Tval),
            false,
            false
        )
    end
    integrator.state.stop = false
    return integrator.state
end

function _forces_and_energy!(system::System, integrator::ConjugateGradient)
    _update_forces!(system)
    return system.forces, integrator.energy(system)
end

function step!(integrator::ConjugateGradient, system::System)
    state = _ensure_state!(integrator, system)
    state.stop = false

    forces, energy = _forces_and_energy!(system, integrator)
    if !state.initialized
        _neg!(state.grad, forces)
        _copy!(state.dir, forces)
        state.energy = energy
        state.initialized = true
    else
        state.energy = energy
    end

    grad = state.grad
    dir = state.dir
    g_norm2 = _vecdot(grad, grad)
    if sqrt(g_norm2) <= integrator.tol
        state.initialized = false
        state.stop = true
        return nothing
    end

    slope = _vecdot(grad, dir)
    if slope ≥ 0
        _copy!(dir, system.forces)
        slope = _vecdot(grad, dir)
        if slope ≥ 0
            state.initialized = false
            state.stop = true
            return nothing
        end
    end

    energy0 = state.energy
    _copy!(state.backup, system.positions)

    α = integrator.alpha0
    accepted = false
    while α > integrator.min_alpha
        @inbounds for i in eachindex(system.positions)
            system.positions[i] = state.backup[i] + α * dir[i]
        end
        integrator.wrap && wrap_positions!(system.positions, system.box)

        forces, trial_energy = _forces_and_energy!(system, integrator)
        if trial_energy <= energy0 + integrator.c1 * α * slope
            state.energy = trial_energy
            accepted = true
            break
        else
            α *= 0.5
        end
    end

    if !accepted
        @inbounds for i in eachindex(system.positions)
            system.positions[i] = state.backup[i]
        end
        integrator.wrap && wrap_positions!(system.positions, system.box)
        _forces_and_energy!(system, integrator)
        state.initialized = false
        state.stop = true
        return nothing
    end

    _neg!(state.grad_tmp, system.forces)
    g_norm2_new = _vecdot(state.grad_tmp, state.grad_tmp)
    if sqrt(g_norm2_new) <= integrator.tol
        state.initialized = false
        state.stop = true
        return nothing
    end

    denom = max(_vecdot(grad, grad), eps(eltype(grad[1])))
    beta_num = zero(eltype(state.grad_tmp[1]))
    @inbounds for i in eachindex(state.grad_tmp)
        diff = state.grad_tmp[i] - grad[i]
        beta_num += dot(state.grad_tmp[i], diff)
    end
    beta = max(beta_num / denom, zero(beta_num))
    _copy!(grad, state.grad_tmp)
    @inbounds for i in eachindex(dir)
        dir[i] = system.forces[i] + beta * dir[i]
    end
    state.stop = false
    return nothing
end

stop_requested(integrator::ConjugateGradient) = (integrator.state isa CGState) && integrator.state.stop

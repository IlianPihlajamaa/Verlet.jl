"""
    ConjugateGradient(energy; tol=1e-8, alpha0=1.0, min_alpha=1e-8, c1=1e-4, wrap=true)

Polak–Ribière non-linear conjugate-gradient minimiser. Forces are obtained from
the system's `ForceField`; the user supplies an `energy(system)` callable used
for the line search. The integrator interprets `nsteps` as the maximum number of
iterations. Set `wrap=true` to apply periodic wrapping after each position
update.
"""
struct ConjugateGradient{T,E} <: AbstractIntegrator
    tol::T
    alpha0::T
    min_alpha::T
    c1::T
    wrap::Bool
    energy::E
end

function ConjugateGradient(energy::Function; tol::Real=1e-8, alpha0::Real=1.0,
                           min_alpha::Real=1e-8, c1::Real=1e-4, wrap::Bool=true)
    T = promote_type(typeof(tol), typeof(alpha0), typeof(min_alpha), typeof(c1))
    ConjugateGradient{T,typeof(energy)}(T(tol), T(alpha0), T(min_alpha), T(c1), wrap, energy)
end

@inline function _vecdot(a, b)
    @assert !isempty(a) "Cannot compute conjugate-gradient direction for empty system"
    s = zero(eltype(a[1]))
    @inbounds for i in eachindex(a)
        s += dot(a[i], b[i])
    end
    return s
end

@inline function _neg!(dest, src)
    @inbounds for i in eachindex(src)
        dest[i] = -src[i]
    end
    return dest
end

@inline function _copy!(dest, src)
    @inbounds for i in eachindex(src)
        dest[i] = src[i]
    end
    return dest
end

function _forces_and_energy!(system::System, integrator::ConjugateGradient)
    _update_forces!(system)
    energy = integrator.energy(system)
    return system.forces, energy
end

function integrate!(integrator::ConjugateGradient, system::System, nsteps::Integer, args...; callback=nothing)
    isempty(args) || throw(ArgumentError("ConjugateGradient does not accept extra positional arguments"))
    nsteps < 0 && throw(ArgumentError("nsteps must be non-negative"))
    nsteps == 0 && return system
    length(system.positions) == 0 && return system

    forces, energy = _forces_and_energy!(system, integrator)
    dir = Vector{typeof(system.positions[1])}(undef, length(forces))
    grad = Vector{typeof(system.positions[1])}(undef, length(forces))
    grad_new = similar(grad)
    _neg!(grad, forces)
    _copy!(dir, forces)
    g_norm2 = _vecdot(grad, grad)
    if sqrt(g_norm2) <= integrator.tol
        if callback !== nothing
            callback(system, 0, energy)
        end
        return system
    end

    positions_backup = copy(system.positions)
    forces_trial = Vector{typeof(system.positions[1])}(undef, length(forces))

    for iter in 1:nsteps
        slope = _vecdot(grad, dir)
        if slope ≥ 0
            _copy!(dir, forces)
            slope = _vecdot(grad, dir)
            slope ≥ 0 && break
        end

        energy0 = energy
        _copy!(positions_backup, system.positions)

        α = integrator.alpha0
        accepted = false
        while α > integrator.min_alpha
            @inbounds for i in eachindex(system.positions)
                system.positions[i] = positions_backup[i] + α * dir[i]
            end
            integrator.wrap && wrap_positions!(system.positions, system.box)

            forces_trial, energy_trial = _forces_and_energy!(system, integrator)
            if energy_trial <= energy0 + integrator.c1 * α * slope
                accepted = true
                energy = energy_trial
                _copy!(forces, forces_trial)
                break
            else
                α *= 0.5
            end
        end

        if !accepted
            @inbounds for i in eachindex(system.positions)
                system.positions[i] = positions_backup[i]
            end
            integrator.wrap && wrap_positions!(system.positions, system.box)
            _forces_and_energy!(system, integrator)
            break
        end

        _neg!(grad_new, forces)
        g_norm2 = _vecdot(grad_new, grad_new)
        if callback !== nothing && callback(system, iter, energy) === false
            break
        end
        if sqrt(g_norm2) <= integrator.tol
            break
        end

        denom = max(_vecdot(grad, grad), eps(eltype(grad[1])))
        beta_num = zero(eltype(grad_new[1]))
        @inbounds for i in eachindex(grad_new)
            diff = grad_new[i] - grad[i]
            beta_num += dot(grad_new[i], diff)
        end
        beta = max(beta_num / denom, zero(beta_num))
        _copy!(grad, grad_new)
        @inbounds for i in eachindex(dir)
            dir[i] = forces[i] + beta * dir[i]
        end
    end

    return system
end

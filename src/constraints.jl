"""
    struct DistanceConstraints

Immutable set of pairwise distance constraints for SHAKE/RATTLE.
"""
struct DistanceConstraints
    i::Vector{Int}
    j::Vector{Int}
    r0::Vector{Float64}
    tol::Float64
    maxiter::Int
    use_minimum_image::Bool
end

function DistanceConstraints(pairs::Vector{<:Tuple}, lengths::Vector{<:Real}; tol=1e-8, maxiter=50, use_minimum_image=true)
    i = Int[i for (i, _) in pairs]
    j = Int[j for (_, j) in pairs]
    r0 = Float64[length for length in lengths]
    return DistanceConstraints(i, j, r0, tol, maxiter, use_minimum_image)
end


function _displacement!(Δ, ri, rj, cons::DistanceConstraints, box=nothing)
    Δ .= ri .- rj
    if cons.use_minimum_image && box !== nothing
        minimum_image!(Δ, box)
    end
    return Δ
end

"""
    apply_shake!(ps, cons, dt)

Iteratively correct positions to satisfy constraints.
"""
function apply_shake!(ps, cons::DistanceConstraints, dt)
    N, D = size(ps.positions)
    Δ = zeros(D)
    masses = ps.masses
    for iter in 1:cons.maxiter
        maxviol = 0.0
        for (k, (i, j)) in enumerate(zip(cons.i, cons.j))
            ri = view(ps.positions, i, :)
            rj = view(ps.positions, j, :)
            _displacement!(Δ, ri, rj, cons)
            dist2 = dot(Δ, Δ)
            C = dist2 - cons.r0[k]^2
            maxviol = max(maxviol, abs(C))
            if abs(C) > cons.tol
                σ = (1/masses[i] + 1/masses[j]) * dist2
                if σ < eps()
                    error("Constraint ill-conditioned: σ≈0 between atoms $i and $j")
                end
                Δλ = -C / (2σ)
                @. ri += (Δλ / masses[i]) * Δ
                @. rj -= (Δλ / masses[j]) * Δ
            end
        end
        if maxviol <= cons.tol
            return
        end
    end
    error("SHAKE did not converge within $(cons.maxiter) iterations")
end

"""
    apply_rattle!(ps, cons)

Correct velocities to satisfy velocity constraints.
"""
function apply_rattle!(ps, cons::DistanceConstraints)
    N, D = size(ps.positions)
    Δ = zeros(D)
    masses = ps.masses
    for iter in 1:cons.maxiter
        maxviol = 0.0
        for (k, (i, j)) in enumerate(zip(cons.i, cons.j))
            ri = view(ps.positions, i, :)
            rj = view(ps.positions, j, :)
            _displacement!(Δ, ri, rj, cons)
            vi = view(ps.velocities, i, :)
            vj = view(ps.velocities, j, :)
            vrel = vi .- vj
            dotdv = dot(Δ, vrel)
            maxviol = max(maxviol, abs(dotdv))
            if abs(dotdv) > cons.tol
                τ = (1/masses[i] + 1/masses[j]) * dot(Δ, Δ)
                μ = -dotdv / τ
                @. vi += (μ / masses[i]) * Δ
                @. vj -= (μ / masses[j]) * Δ
            end
        end
        if maxviol <= cons.tol
            return
        end
    end
    error("RATTLE did not converge within $(cons.maxiter) iterations")
end

"""
    velocity_verlet_shake_rattle!(ps, forces, dt, cons)

Constrained velocity Verlet step with SHAKE/RATTLE.
"""
function velocity_verlet_shake_rattle!(ps, forces, dt, cons::DistanceConstraints)
    invm = 1.0 ./ ps.masses
    # Half kick
    F = forces(ps.positions)
    ps.velocities .+= 0.5 .* dt .* F .* invm
    # Drift
    ps.positions .+= dt .* ps.velocities
    # SHAKE
    apply_shake!(ps, cons, dt)
    # New forces
    F = forces(ps.positions)
    # Half kick
    ps.velocities .+= 0.5 .* dt .* F .* invm
    # RATTLE
    apply_rattle!(ps, cons)
    return ps
end

"""
    remove_com_motion!(ps; which=:velocity)

Remove center-of-mass motion from velocities/positions.
"""
function remove_com_motion!(ps; which=:velocity)
    m = ps.masses
    Mtot = sum(m)
    if which == :velocity || which == :both
        Vcom = (ps.velocities' * m) / Mtot
        ps.velocities .-= Vcom'
    end
    if which == :position || which == :both
        Rcom = (ps.positions' * m) / Mtot
        ps.positions .-= Rcom'
    end
    return ps
end

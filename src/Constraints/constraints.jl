"""
    struct DistanceConstraints{T_int,T_float}

Immutable set of pairwise distance constraints for SHAKE/RATTLE.
"""
struct DistanceConstraints{T_int,T_float}
    i::Vector{T_int}
    j::Vector{T_int}
    r0::Vector{T_float}
    tol::T_float
    maxiter::T_int
    use_minimum_image::Bool
end

function DistanceConstraints(pairs::Vector{<:Tuple}, lengths::Vector{<:Real}; tol=1e-8, maxiter=50, use_minimum_image=true)
    i = T_int[i for (i, _) in pairs]
    j = T_int[j for (_, j) in pairs]
    r0 = T_float[length for length in lengths]
    return DistanceConstraints{T_int,T_float}(i, j, r0, T_float(tol), T_int(maxiter), use_minimum_image)
end


function _displacement!(ri, rj, cons::DistanceConstraints, box=nothing)
    Δnew = ri - rj
    if cons.use_minimum_image && box !== nothing
        Δnew = minimum_image(Δnew, box)
    end
    return Δnew
end

"""
    apply_shake!(ps, cons, dt)

Iteratively correct positions to satisfy constraints.
"""
function apply_shake!(ps, cons::DistanceConstraints, dt)
    masses = ps.masses
    for iter in 1:cons.maxiter
        maxviol = 0.0
        for (k, (i, j)) in enumerate(zip(cons.i, cons.j))
            ri = ps.positions[i]
            rj = ps.positions[j]
            Δ = _displacement!(ri, rj, cons)
            dist2 = dot(Δ, Δ)
            C = dist2 - cons.r0[k]^2
            maxviol = max(maxviol, abs(C))
            if abs(C) > cons.tol
                σ = (1/masses[i] + 1/masses[j]) * dist2
                if σ < eps()
                    error("Constraint ill-conditioned: σ≈0 between atoms $i and $j")
                end
                Δλ = -C / (2σ)
                ps.positions[i] += (Δλ / masses[i]) * Δ
                ps.positions[j] -= (Δλ / masses[j]) * Δ
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
    masses = ps.masses
    for iter in 1:cons.maxiter
        maxviol = 0.0
        for (k, (i, j)) in enumerate(zip(cons.i, cons.j))
            ri = ps.positions[i]
            rj = ps.positions[j]
            Δ = _displacement!(ri, rj, cons)
            vi = ps.velocities[i]
            vj = ps.velocities[j]
            vrel = vi - vj
            dotdv = dot(Δ, vrel)
            maxviol = max(maxviol, abs(dotdv))
            if abs(dotdv) > cons.tol
                τ = (1/masses[i] + 1/masses[j]) * dot(Δ, Δ)
                μ = -dotdv / τ
                ps.velocities[i] += (μ / masses[i]) * Δ
                ps.velocities[j] -= (μ / masses[j]) * Δ
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
    ps.velocities .= ps.velocities .+ 0.5 .* dt .* F .* invm
    # Drift
    ps.positions .= ps.positions .+ dt .* ps.velocities
    # SHAKE
    apply_shake!(ps, cons, dt)
    # New forces
    F = forces(ps.positions)
    # Half kick
    ps.velocities .= ps.velocities .+ 0.5 .* dt .* F .* invm
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
        Vcom = sum(ps.velocities[i] * m[i] for i in eachindex(m)) / Mtot
        ps.velocities .= ps.velocities .- Ref(Vcom)
    end
    if which == :position || which == :both
        Rcom = sum(ps.positions[i] * m[i] for i in eachindex(m)) / Mtot
        ps.positions .= ps.positions .- Ref(Rcom)
    end
    return ps
end

# ------------------------------------------------------------
# Constraint residual monitoring (docstring and example)
# ------------------------------------------------------------
"""
    constraint_residuals(ps::ParticleSystem, cons::DistanceConstraints) -> (; maxC, rmsC, maxCd, rmsCd)

Compute the constraint residuals for a system subject to pairwise distance constraints.

For each constraint `l` with atoms `(i, j)` and target distance `r0_l`:

- Position residual: `C_l = ||r_i - r_j||^2 - r0_l^2`
- Velocity residual: `Ċ_l = 2 (r_i - r_j) ⋅ (v_i - v_j)`

Returns a named tuple with `maxC`, `rmsC`, `maxCd`, `rmsCd`.

Example
```julia
ps = ParticleSystem([SVector(0.0, 0, 0), SVector(1.0, 0, 0)], [SVector(0.0, 0, 0), SVector(0.0, 0, 0)], [1.0, 1.0])
cons = DistanceConstraints([(1,2)], [1.0])
constraint_residuals(ps, cons)
```
"""
function constraint_residuals(ps, cons::DistanceConstraints)
    R = ps.positions
    V = ps.velocities

    # Use the actual fields of DistanceConstraints
    is = cons.i
    js = cons.j
    lengths = cons.r0

    L = length(is)
    maxC  = 0.0
    maxCd = 0.0
    accC2  = 0.0
    accCd2 = 0.0

    for ℓ in 1:L
        i = is[ℓ]
        j = js[ℓ]
        ri = R[i]
        rj = R[j]
        d = _displacement!(ri, rj, cons)
        vij = V[i] - V[j]

        r0 = lengths[ℓ]
        Cl = dot(d, d) - r0^2
        Cdl = 2.0 * dot(d, vij)

        maxC  = max(maxC,  abs(Cl))
        maxCd = max(maxCd, abs(Cdl))
        accC2  += Cl^2
        accCd2 += Cdl^2
    end

    rmsC  = L == 0 ? 0.0 : sqrt(accC2  / L)
    rmsCd = L == 0 ? 0.0 : sqrt(accCd2 / L)
    return (; maxC, rmsC, maxCd, rmsCd)
end


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
    KE = T(0.5) * sum(m[i] * sum(abs2.(v[i])) for i in eachindex(v))
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
    sys.velocities .*= λ
    return sys
end

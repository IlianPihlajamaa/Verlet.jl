function _ensure_forcefield(system::System)
    ff = system.forcefield
    ff === nothing && throw(ArgumentError("system has no forcefield assigned"))
    return ff
end

function _update_forces!(system::System)
    ff = _ensure_forcefield(system)
    compute_all_forces!(system, ff)
    return system.forces
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

"""
    struct CubicBox{T<:Real}

Simple cubic periodic box with side length `L`.
"""
struct CubicBox{T_float}
    L::T_float
end

"""
    minimum_image!(Δ, box::CubicBox)

Apply the minimum-image convention to the displacement vector `Δ` in-place
for the periodic `box` (wraps components to (-L/2, L/2]).
"""
function minimum_image!(Δ::AbstractVector, box::CubicBox)
    half = box.L / 2
    @inbounds for k in eachindex(Δ)
        if Δ[k] >  half
            Δ[k] -= box.L
        elseif Δ[k] <= -half
            Δ[k] += box.L
        end
    end
    return Δ
end

"""
    wrap_positions!(R, box::CubicBox)

Wrap particle positions `R` into the primary periodic image of `box` in-place.
This is useful before building neighbor lists or measuring displacements.
"""
function wrap_positions!(R::AbstractMatrix, box::CubicBox)
    half = box.L / 2
    L = box.L
    @inbounds for i in axes(R,1), k in axes(R,2)
        while R[i,k] >  half
            R[i,k] -= L
        end
        while R[i,k] <= -half
            R[i,k] += L
        end
    end
    return R
end


function box_length(box::CubicBox)
    return box.L
end
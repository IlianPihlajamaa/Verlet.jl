"""
    abstract type AbstractBox{T<:AbstractFloat} end

Abstract type for simulation boxes.
"""
abstract type AbstractBox{T<:AbstractFloat} end

"""
    struct CubicBox{T<:AbstractFloat} <: AbstractBox{T}

Simple cubic periodic box with side length `L`.
"""
struct CubicBox{T<:AbstractFloat} <: AbstractBox{T}
    L::T
end


"""
    minimum_image(Δ, box::CubicBox)
    minimum_image(Δ, L)

Apply the minimum-image convention to the displacement vector `Δ`
for the periodic `box` (wraps components to (-L/2, L/2]) and return the new vector.
"""
@inline minimum_image(Δ::AbstractVector, box::CubicBox) = minimum_image(Δ::AbstractVector, box_length(box))

@inline function minimum_image(Δ::AbstractVector, L)
    Δnew = Δ - L * round.(Δ / L)
    return Δnew
end



"""
    wrap_positions!(R, box::CubicBox)

Wrap particle positions `R` into the primary periodic image of `box` in-place.
This is useful before building neighbor lists or measuring displacements.

The resulting positions will be in the range (-L/2, L/2].
"""
function wrap_positions!(R::Vector{SVector{Dims, T}}, box::CubicBox{T}) where {Dims,T<:AbstractFloat}
    for i in eachindex(R)
        R[i] = minimum_image(R[i], box)
    end

    return R
end


@inline function box_length(box::CubicBox)
    return box.L
end
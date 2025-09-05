"""
    struct CubicBox{T<:Real}

Simple cubic periodic box with side length `L`.
"""
struct CubicBox{T<:Real}
    L::T
end

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

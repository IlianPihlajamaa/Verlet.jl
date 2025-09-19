maybe_detach(value) = Base.ismutabletype(typeof(value)) ? _copy_or_deepcopy(value) : value

function _copy_or_deepcopy(value)
    if hasmethod(copy, Tuple{typeof(value)})
        return copy(value)
    else
        return deepcopy(value)
    end
end

function materialize_value(::Type{T}, value) where {T}
    converted = convert(T, value)
    return maybe_detach(converted)
end

function default_bilinear(::Type{TA}, ::Type{TB}) where {TA,TB}
    if TA <: Number && TB <: Number
        return (a::TA, b::TB) -> a * b
    elseif hasmethod(dot, Tuple{TA,TB})
        return (a::TA, b::TB) -> dot(a, b)
    else
        throw(ArgumentError("Provide `product` keyword argument to specify how to contract $(TA) and $(TB)`"))
    end
end

mutable struct SumStorage{T}
    value::Union{Nothing,T}
end

SumStorage{T}() where {T} = SumStorage{T}(nothing)

_copy_value(x::AbstractArray) = copy(x)
_copy_value(x) = x

function _add_values(acc::AbstractArray, value)
    if Base.ismutabletype(typeof(acc))
        @inbounds @. acc += value
        return acc
    else
        return acc .+ value
    end
end

function _add_values(acc, value)
    return acc + value
end

function _scale_value(value::AbstractArray, factor::Int)
    return value ./ factor
end

function _scale_value(value, factor::Int)
    return value / factor
end

function _persist_for_sum(::Type{T}, value) where {T}
    return materialize_value(T, value)
end

function _accumulate!(storage::SumStorage{T}, value) where {T}
    if storage.value === nothing
        storage.value = _persist_for_sum(T, value)
    else
        current = storage.value::T
        storage.value = _add_values(current, convert(T, value))
    end
    return storage
end

function _mean_from_sum(storage::SumStorage, nsamples::Int)
    storage.value === nothing && return nothing
    return _scale_value(_copy_value(storage.value), nsamples)
end

mutable struct StaticCorrelationLogger{ObsTuple<:Tuple{Vararg{Observable}}, TypeTuple<:Tuple, SumTuple<:Tuple, MatrixT<:AbstractMatrix, F} <: AbstractLogger
    observables::ObsTuple
    schedule::LoggingSchedule
    nsamples::Int
    value_types::TypeTuple
    sums::SumTuple
    sum_products::MatrixT
    product::F
end

function StaticCorrelationLogger(observables::Tuple{Vararg{Observable}}; when_to_log::AbstractVector{<:Integer}=Int[], stride::Integer=1, product=nothing)
    schedule = LoggingSchedule(; stride, when_to_log)
    n = length(observables)
    value_types = ntuple(i -> out_type(observables[i]), n)
    sums = ntuple(i -> SumStorage{value_types[i]}(), n)
    TA = value_types[1]
    TB = n >= 2 ? value_types[2] : TA
    prod_fn = product === nothing ? default_bilinear(TA, TB) : product
    TC = Base.promote_op(prod_fn, TA, TB)
    sum_products = fill(zero(TC), n, n)
    return StaticCorrelationLogger(observables, schedule, 0, value_types, sums, sum_products, prod_fn)
end

StaticCorrelationLogger(observable::Observable; kwargs...) = StaticCorrelationLogger((observable,); kwargs...)

when_to_log(logger::StaticCorrelationLogger) = when_to_log(logger.schedule)
stride(logger::StaticCorrelationLogger) = stride(logger.schedule)

function step!(logger::StaticCorrelationLogger, system::System, step::Integer)
    if !should_sample!(logger.schedule, step)
        return logger
    end

    N = length(logger.observables)
    values = ntuple(i -> _persist_for_sum(logger.value_types[i], observe!(logger.observables[i], system, step)), Val(N))

    sums = logger.sums
    for i in 1:N
        _accumulate!(sums[i], values[i])
    end

    for i in 1:N
        vi = values[i]
        for j in i:N
            inner = logger.product(vi, values[j])
            logger.sum_products[i, j] += inner
            if i != j
                logger.sum_products[j, i] += inner
            end
        end
    end

    logger.nsamples += 1
    return logger
end

function means(logger::StaticCorrelationLogger)
    N = length(logger.sums)
    logger.nsamples == 0 && return ntuple(_ -> nothing, Val(N))
    return ntuple(i -> _mean_from_sum(logger.sums[i], logger.nsamples), Val(N))
end

function covars(logger::StaticCorrelationLogger)
    logger.nsamples == 0 && return fill(zero(eltype(logger.sum_products)), size(logger.sum_products))
    return logger.sum_products ./ logger.nsamples
end

function data(logger::StaticCorrelationLogger)
    return (; nsamples = logger.nsamples, means = means(logger), covars = covars(logger))
end

function save(io::IO, logger::StaticCorrelationLogger)
    println(io, "# nsamples\t", logger.nsamples)
    println(io, "# means")
    μ = means(logger)
    for (obs, m) in zip(logger.observables, μ)
        println(io, string(observed_quantity(obs), '\t', repr(m)))
    end
    println(io, "# covariances")
    cov = covars(logger)
    for i in 1:size(cov, 1)
        println(io, join(cov[i, :], '\t'))
    end
    return nothing
end

function save(path::AbstractString, logger::StaticCorrelationLogger)
    open(path, "w") do io
        save(io, logger)
    end
    return path
end

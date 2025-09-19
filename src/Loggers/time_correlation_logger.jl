mutable struct WindowBuffer{TA,TB,TC}
    stride::Int
    length::Int
    filled::Int
    write_index::Int
    a_values::Vector{TA}
    b_values::Vector{TB}
    corr_sum::Vector{TC}
    counts::Vector{Int}
end

function WindowBuffer{TA,TB,TC}(stride::Int, length::Int) where {TA,TB,TC}
    a_values = Vector{TA}(undef, length)
    b_values = Vector{TB}(undef, length)
    corr_sum = Vector{TC}(undef, length)
    fill!(corr_sum, zero(TC))
    counts = zeros(Int, length)
    return WindowBuffer{TA,TB,TC}(stride, length, 0, 0, a_values, b_values, corr_sum, counts)
end

function _update_window!(window::WindowBuffer{TA,TB,TC}, sampleA::TA, sampleB::TB, product) where {TA,TB,TC}
    len = window.length
    new_index = window.write_index == len ? 1 : window.write_index + 1
    window.a_values[new_index] = sampleA
    window.b_values[new_index] = sampleB
    if window.filled < len
        window.filled += 1
    end
    window.write_index = new_index

    origin = new_index - window.filled + 1
    if origin <= 0
        origin += len
    end

    aval = window.a_values[origin]
    idx_target = origin
    for lag in 0:window.filled-1
        if idx_target > len
            idx_target -= len
        end
        bval = window.b_values[idx_target]
        window.corr_sum[lag + 1] += product(aval, bval)
        window.counts[lag + 1] += 1
        idx_target += 1
    end

    return nothing
end

mutable struct TimeCorrelationLogger{ObsTuple<:Tuple{Vararg{Observable}},TA,TB,TC,F} <: AbstractLogger
    observables::ObsTuple
    block_length::Int
    blocks::Vector{WindowBuffer{TA,TB,TC}}
    timestep::Float64
    sample_stride::Int
    product::F
end

function _check_out_type(name::Symbol, T::Type)
    T === Any && throw(ArgumentError("observable $(name) must declare a concrete out_type"))
    return T
end

function TimeCorrelationLogger(observables::Tuple{Vararg{Observable}}; block_length::Integer=32, nblocks::Integer=6, min_stride::Integer=1, timestep::Real=1.0, product=nothing)
    length(observables) in (1, 2) || throw(ArgumentError("TimeCorrelationLogger expects one (auto) or two (cross) observables"))
    min_stride < 1 && throw(ArgumentError("min_stride must be ≥ 1"))
    block_length < 2 && throw(ArgumentError("block_length must be ≥ 2"))
    nblocks < 1 && throw(ArgumentError("nblocks must be ≥ 1"))

    obsA = observables[1]
    TA = _check_out_type(observed_quantity(obsA), out_type(obsA))
    if length(observables) == 2
        obsB = observables[2]
        TB = _check_out_type(observed_quantity(obsB), out_type(obsB))
    else
        TB = TA
    end

    prod_fn = product === nothing ? default_bilinear(TA, TB) : product
    TC = Base.promote_op(prod_fn, TA, TB)

    blocks = Vector{WindowBuffer{TA,TB,TC}}(undef, nblocks)
    for b in 1:nblocks
        stride = min_stride * block_length^(b - 1)
        blocks[b] = WindowBuffer{TA,TB,TC}(stride, block_length)
    end

    return TimeCorrelationLogger{typeof(observables),TA,TB,TC,typeof(prod_fn)}(observables, Int(block_length), blocks, Float64(timestep), Int(min_stride), prod_fn)
end

TimeCorrelationLogger(obs::Observable; kwargs...) = TimeCorrelationLogger((obs,); kwargs...)
TimeCorrelationLogger(obsA::Observable, obsB::Observable; kwargs...) = TimeCorrelationLogger((obsA, obsB); kwargs...)

function _collect_sample(logger::TimeCorrelationLogger{ObsTuple,TA,TB,TC,F}, system::System, step::Integer) where {ObsTuple,TA,TB,TC,F}
    obsA = logger.observables[1]
    rawA = observe!(obsA, system, step)
    sampleA = materialize_value(TA, rawA)

    sampleB = if length(logger.observables) == 1
        sampleA
    else
        rawB = observe!(logger.observables[2], system, step)
        materialize_value(TB, rawB)
    end

    return sampleA, sampleB
end

function step!(logger::TimeCorrelationLogger, system::System, step::Integer)
    if step % logger.sample_stride != 0
        return logger
    end

    sampleA, sampleB = _collect_sample(logger, system, step)

    for block in logger.blocks
        if step % block.stride == 0
            _update_window!(block, sampleA, sampleB, logger.product)
        end
    end

    return logger
end

function _collect_block_data(block::WindowBuffer, timestep::Float64)
    filled = block.filled
    if filled == 0
        return Float64[], similar(block.corr_sum, 0), Int[]
    end
    stride = block.stride
    times = Vector{Float64}(undef, filled)
    values = Vector{eltype(block.corr_sum)}(undef, filled)
    counts = Vector{Int}(undef, filled)
    for lag in 0:(filled - 1)
        times[lag + 1] = stride * lag * timestep
        count = block.counts[lag + 1]
        counts[lag + 1] = count
        values[lag + 1] = count > 0 ? block.corr_sum[lag + 1] / count : zero(eltype(block.corr_sum))
    end
    return times, values, counts
end

function data(logger::TimeCorrelationLogger)
    all_times = Float64[]
    value_type = eltype(first(logger.blocks).corr_sum)
    all_values = Vector{value_type}()
    all_counts = Int[]
    strides = Int[]

    for (idx, block) in enumerate(logger.blocks)
        times, values, counts = _collect_block_data(block, logger.timestep)
        if idx > 1 && length(times) > 1
            times = times[2:end]
            values = values[2:end]
            counts = counts[2:end]
        elseif idx > 1
            times = Float64[]
            values = Float64[]
            counts = Int[]
        end
        append!(all_times, times)
        append!(all_values, values)
        append!(all_counts, counts)
        push!(strides, block.stride)
    end

    return (; times = all_times, values = all_values, counts = all_counts, block_strides = strides, block_length = logger.block_length, sample_stride = logger.sample_stride)
end

function save(io::IO, logger::TimeCorrelationLogger)
    info = data(logger)
    println(io, "# sample_stride: ", info.sample_stride)
    println(io, "# block length: ", logger.block_length)
    println(io, "# block strides: ", join(info.block_strides, ", "))
    println(io, "time\tvalue\tcount")
    for (t, v, c) in zip(info.times, info.values, info.counts)
        println(io, string(t, '\t', repr(v), '\t', c))
    end
    return nothing
end

function save(path::AbstractString, logger::TimeCorrelationLogger)
    open(path, "w") do io
        save(io, logger)
    end
    return path
end

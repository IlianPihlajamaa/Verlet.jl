abstract type AbstractLogger end

mutable struct ObservableLogger{ObsT<:Observable,T} <: AbstractLogger
    observable::ObsT
    schedule::LoggingSchedule
    steps::Vector{Int}
    values::Vector{T}
    max_samples::Union{Nothing,Int}
end

function ObservableLogger(obs::Observable; when_to_log::AbstractVector{<:Integer}=Int[], stride::Integer=1, max_samples::Union{Nothing,Integer}=nothing)
    schedule = LoggingSchedule(; stride, when_to_log)
    T = out_type(obs)
    value_buffer = Vector{T}()
    max_samples !== nothing && max_samples < 1 && throw(ArgumentError("max_samples must be ≥ 1"))
    return ObservableLogger(obs, schedule, Int[], value_buffer, max_samples === nothing ? nothing : Int(max_samples))
end

when_to_log(logger::ObservableLogger) = when_to_log(logger.schedule)
stride(logger::ObservableLogger) = stride(logger.schedule)

function Base.length(logger::ObservableLogger)
    return length(logger.steps)
end

function step!(logger::ObservableLogger{ObsT,T}, system::System, step::Integer) where {ObsT<:Observable,T}
    if should_sample!(logger.schedule, step)
        sample = observe!(logger.observable, system, step)
        stored = materialize_value(T, sample)
        push!(logger.steps, step)
        push!(logger.values, stored)
        if logger.max_samples !== nothing && length(logger.steps) > logger.max_samples
            popfirst!(logger.steps)
            popfirst!(logger.values)
        end
    end
    return logger
end

function data(logger::ObservableLogger)
    return (; steps = copy(logger.steps), values = copy(logger.values))
end

function save(io::IO, logger::ObservableLogger)
    for (step, value) in zip(logger.steps, logger.values)
        print(io, step)
        print(io, '\t')
        if value isa AbstractArray
            println(io, join(value, '\t'))
        else
            println(io, repr(value))
        end
    end
    return nothing
end

function save(path::AbstractString, logger::ObservableLogger)
    open(path, "w") do io
        save(io, logger)
    end
    return path
end

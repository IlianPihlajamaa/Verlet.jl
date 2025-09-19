mutable struct LoggingSchedule
    stride::Int
    when_to_log::Vector{Int}
    next_index::Int
end

function LoggingSchedule(; stride::Integer=1, when_to_log::AbstractVector{<:Integer}=Int[])
    stride < 0 && throw(ArgumentError("stride must be non-negative"))
    steps = sort!(unique!(collect(Int, when_to_log)))
    return LoggingSchedule(Int(stride), steps, 1)
end

function Base.copy(schedule::LoggingSchedule)
    LoggingSchedule(schedule.stride, copy(schedule.when_to_log), schedule.next_index)
end

function should_sample!(schedule::LoggingSchedule, step::Integer)
    stride_hit = schedule.stride > 0 && step % schedule.stride == 0

    when_len = length(schedule.when_to_log)
    while schedule.next_index <= when_len && schedule.when_to_log[schedule.next_index] < step
        schedule.next_index += 1
    end

    explicit_hit = schedule.next_index <= when_len && schedule.when_to_log[schedule.next_index] == step
    if explicit_hit
        schedule.next_index += 1
    end

    return stride_hit || explicit_hit
end

when_to_log(schedule::LoggingSchedule) = schedule.when_to_log
stride(schedule::LoggingSchedule) = schedule.stride

using Test
using Verlet
import Verlet: observe!, out_type, observed_quantity
using StaticArrays

function sample_system()
    positions = [SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(1.0, 0.0, 0.0)]
    velocities = [SVector{3,Float64}(0.1, 0.0, 0.0), SVector{3,Float64}(0.0, 0.2, 0.0)]
    forces = [SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0)]
    masses = [1.0, 2.0]
    box = CubicBox(10.0)
    types = ones(Int, 2)
    type_names = Dict(1 => :A)
    ff = ForceField(())
    return System(positions, velocities, forces, masses, box, types, type_names; forcefield=ff)
end

@testset "ObservableLogger basics" begin
    sys = sample_system()
    obs = DensityObservable()
    logger = ObservableLogger(obs; stride=2, when_to_log=[3, 4], max_samples=3)

    for step in 1:6
        step!(logger, sys, step)
    end

    logged = data(logger)
    @test logged.steps == [3, 4, 6]
    expected_density = Float64(natoms(sys)) / volume(sys)
    @test all(isapprox(val, expected_density; rtol=1e-12) for val in logged.values)
    @test length(logged.steps) == 3
end

@testset "Default observables" begin
    sys = sample_system()
    @test observe!(TemperatureObservable(), sys, 1) > 0
    @test observe!(KineticEnergyObservable(), sys, 1) ≈ kinetic_energy(sys)
    @test observe!(PotentialEnergyObservable(), sys, 1) ≈ 0.0
    @test observe!(TotalEnergyObservable(), sys, 1) ≈ kinetic_energy(sys)
    velocities = observe!(VelocityObservable(), sys, 1)
    @test velocities isa Vector{Float64}
    forces = observe!(ForceObservable(), sys, 1)
    @test forces isa Vector{Float64}
    @test observe!(VolumeObservable(), sys, 1) ≈ volume(sys)
    pressure = observe!(PressureObservable(), sys, 1)
    @test pressure > 0
end

@testset "StaticCorrelationLogger" begin
    sys = sample_system()
    obs_tuple = (KineticEnergyObservable(), DensityObservable())
    logger = StaticCorrelationLogger(obs_tuple)

    for step in 1:5
        step!(logger, sys, step)
    end

    @test logger.nsamples == 5
    μ = means(logger)
    @test length(μ) == 2
    @test μ[1] ≈ kinetic_energy(sys)
    @test μ[2] ≈ Float64(natoms(sys)) / volume(sys)
    cov = covars(logger)
    @test size(cov) == (2, 2)
@test cov[1, 2] ≈ μ[1] * μ[2]

    io = IOBuffer()
    @test save(io, logger) === nothing
end

struct ComplexObservable <: Observable end
observe!(::ComplexObservable, ::System, step::Integer) = ComplexF64(step, -step)
out_type(::ComplexObservable) = ComplexF64
observed_quantity(::ComplexObservable) = :complex_obs

@testset "StaticCorrelationLogger preserves type" begin
    sys = sample_system()
    logger = StaticCorrelationLogger((ComplexObservable(), ComplexObservable()); product = (a, b) -> a * conj(b))
    for step in 1:5
        step!(logger, sys, step)
    end
    cov = covars(logger)
    @test eltype(cov) == ComplexF64
    @test cov[1, 1] isa ComplexF64
end

struct StepObservable <: Observable end
observe!(::StepObservable, ::System, step::Integer) = Float64(step)
out_type(::StepObservable) = Float64
observed_quantity(::StepObservable) = :step

struct SquareObservable <: Observable end
observe!(::SquareObservable, ::System, step::Integer) = Float64(step^2)
out_type(::SquareObservable) = Float64
observed_quantity(::SquareObservable) = :step2

function simulate_multiple_windows(values_a::Vector{Float64}, values_b::Vector{Float64}, block_length::Int, stride::Int, sample_stride::Int)
    len = block_length
    a_buf = Vector{Float64}(undef, len)
    filled = 0
    write_index = 0
    corr_sum = zeros(Float64, len)
    counts = zeros(Int, len)
    for step in 0:(length(values_a) - 1)
        step % sample_stride == 0 || continue
        step % stride == 0 || continue
        write_index = write_index == len ? 1 : write_index + 1
        sample_idx = step + 1
        a_buf[write_index] = values_a[sample_idx]
        if filled < len
            filled += 1
        end
        idx = write_index
        for lag in 0:(filled - 1)
            corr_sum[lag + 1] += a_buf[idx] * values_b[sample_idx]
            counts[lag + 1] += 1
            idx -= 1
            if idx == 0
                idx = len
            end
        end
    end
    return corr_sum, counts, filled
end

@testset "TimeCorrelationLogger (auto)" begin
    sys = sample_system()
    logger = TimeCorrelationLogger(StepObservable(); block_length=3, nblocks=2, min_stride=1, timestep=1.0)
    samples = Float64[]
    for step in 0:12
        push!(samples, Float64(step))
        step!(logger, sys, step)
    end

    for block in logger.blocks
        expected_sum, expected_counts, expected_filled = simulate_multiple_windows(samples, samples, logger.block_length, block.stride, logger.sample_stride)
        @test block.filled == expected_filled
        for lag in 1:block.filled
            @test block.counts[lag] == expected_counts[lag]
            @test block.corr_sum[lag] ≈ expected_sum[lag]
        end
    end

    info = data(logger)
    @test !isempty(info.times)
    @test info.times[1] == 0.0
    @test all(c > 0 for c in info.counts)
end

@testset "TimeCorrelationLogger (cross)" begin
    sys = sample_system()
    logger = TimeCorrelationLogger(StepObservable(), SquareObservable(); block_length=4, nblocks=2, min_stride=2, timestep=0.5)
    samples_a = Float64[]
    samples_b = Float64[]
    for step in 0:20
        push!(samples_a, Float64(step))
        push!(samples_b, Float64(step^2))
        step!(logger, sys, step)
    end

    for block in logger.blocks
        expected_sum, expected_counts, expected_filled = simulate_multiple_windows(samples_a, samples_b, logger.block_length, block.stride, logger.sample_stride)
        @test block.filled == expected_filled
        for lag in 1:block.filled
            @test block.counts[lag] == expected_counts[lag]
            @test block.corr_sum[lag] ≈ expected_sum[lag]
        end
    end

    info = data(logger)
    @test all(info.times .>= 0)
    @test length(info.times) == length(info.values) == length(info.counts)
end

struct ComplexStepObservable <: Observable end
observe!(::ComplexStepObservable, ::System, step::Integer) = ComplexF64(step, step)
out_type(::ComplexStepObservable) = ComplexF64
observed_quantity(::ComplexStepObservable) = :complex_step

@testset "TimeCorrelationLogger preserves value type" begin
    sys = sample_system()
    logger = TimeCorrelationLogger(ComplexStepObservable(); block_length=3, nblocks=2, min_stride=1, timestep=1.0, product=(a, b) -> a * conj(b))
    for step in 0:10
        step!(logger, sys, step)
    end
    info = data(logger)
    @test eltype(info.values) == ComplexF64
    @test info.values[1] isa ComplexF64
end

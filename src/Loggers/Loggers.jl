module Loggers

using LinearAlgebra: dot
using StaticArrays

using ..Core: System, Observable, AbstractBox, natoms, kinetic_energy,
                 compute_potential_energy, volume, system_volume, Dims
import ..Core: observe!, out_type, observed_quantity
using ..Core: AbstractBox
using ..Thermostats: instantaneous_temperature

include("schedule.jl")
include("utils.jl")
include("observable_logger.jl")
include("static_correlation_logger.jl")
include("time_correlation_logger.jl")
include("default_observables.jl")

export AbstractLogger, ObservableLogger, StaticCorrelationLogger
export data, save, step!, when_to_log, stride, means, covars
export DensityObservable, TemperatureObservable, KineticEnergyObservable,
       PotentialEnergyObservable, TotalEnergyObservable, VelocityObservable,
       ForceObservable, VolumeObservable, PressureObservable
export TimeCorrelationLogger

end

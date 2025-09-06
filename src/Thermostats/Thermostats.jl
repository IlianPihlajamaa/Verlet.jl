module Thermostats
using LinearAlgebra, StaticArrays

using ..Core
include("baoab.jl")

# Export thermostat API
export degrees_of_freedom, instantaneous_temperature, velocity_rescale!, langevin_baoab!
export langevin_baoab_constrained!
end

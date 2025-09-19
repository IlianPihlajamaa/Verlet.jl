module Thermostats

using ..Core
using ..Constraints
using LinearAlgebra, StaticArrays, Random

include("baoab.jl")

export degrees_of_freedom, instantaneous_temperature, velocity_rescale!

end

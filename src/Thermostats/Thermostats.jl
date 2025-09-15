module Thermostats

using ..Core
using ..Constraints
using LinearAlgebra, StaticArrays, Random

include("baoab.jl")

export degrees_of_freedom, instantaneous_temperature, velocity_rescale!, langevin_baoab!, langevin_baoab_constrained!

end

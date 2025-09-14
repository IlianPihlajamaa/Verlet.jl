module Thermostats

using ..Core
using ..Constraints
using LinearAlgebra, StaticArrays

include("baoab.jl")

export langevin_baoab!, degrees_of_freedom, instantaneous_temperature, velocity_rescale!, langevin_baoab_constrained!

end

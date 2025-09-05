using Test
using Verlet

using LinearAlgebra

# Core unit tests
include("test_core.jl")

# Neighborlist and cellgrid acceptance tests
include("test_neighborlist.jl")
include("test_cellgrid.jl")

# Thermostat / NVT tests
include("test_thermostats.jl")

# Constraints tests
include("test_constraints.jl")


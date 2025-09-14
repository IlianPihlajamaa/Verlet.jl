import Pkg
Pkg.activate(joinpath(@__DIR__, "."))
Pkg.develop(path=joinpath(@__DIR__, ".."))
Pkg.instantiate()

include(joinpath(@__DIR__, "runtests.jl"))

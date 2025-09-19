# from root: julia --project=docs docs/make.jl

using Documenter
using Verlet

# DocMeta.setdocmeta!(Verlet, :DocTestSetup, :(using Verlet); recursive=true)

makedocs(
    sitename = "Verlet.jl",
    modules = [Verlet, Verlet.Core, Verlet.Neighbors, Verlet.Potentials, Verlet.Constraints, Verlet.Thermostats, Verlet.Loggers],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", "false") == "true"),
    checkdocs = :public,
    pages = [
        "Home" => "index.md",
        "Tutorials" => Any[
            "Tutorial 1 · Hello, Velocity Verlet" => "tutorials/quickstart.md",
            "Tutorial 2 · Pair Potentials & Neighbour Lists" => "tutorials/pair_potentials.md",
            "Tutorial 3 · Constraints in Practice" => "tutorials/constraints.md",
            "Tutorial 4 · Thermostatted Dynamics" => "tutorials/thermostat.md",
            "Tutorial 5 · Observables & Loggers" => "tutorials/observables_pair_distribution.md",
        ],
        "How-to Guides" => Any[
            "Building Systems" => "guide/system.md",
            "Forces & Potentials" => "guide/forces.md",
            "Bonded Interactions" => "guide/bonded.md",
            "Neighbour Lists" => "guide/neighbors.md",
            "Constrained Dynamics" => "guide/constraints.md",
            "Numerical Notes" => "guide/numerics.md",
        ],
        "Reference" => Any[
            "Overview" => "specs/Spec_Module_Index.md",
            "Top-level (Verlet)" => "specs/Spec_Verlet.md",
            "Core" => "specs/Spec_Core.md",
            "Neighbors" => "specs/Spec_Neighbors.md",
            "Potentials" => "specs/Spec_Potentials.md",
            "Constraints" => "specs/Spec_Constraints.md",
            "Thermostats" => "specs/Spec_Thermostats.md",
            "Electrostatics" => "specs/Spec_Electrostatics.md",
            "Loggers" => "specs/Spec_Loggers.md",
        ],
        "API Reference" => "api.md",
    ],
)


deploydocs(; repo = "https://github.com/IlianPihlajamaa/Verlet.jl.git ", devbranch = "main")

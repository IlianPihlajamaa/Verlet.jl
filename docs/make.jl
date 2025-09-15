# from root: julia --project=docs docs/make.jl

using Documenter
using Verlet

# DocMeta.setdocmeta!(Verlet, :DocTestSetup, :(using Verlet); recursive=true)

makedocs(
    sitename = "Verlet.jl",
    modules = [Verlet, Verlet.Core, Verlet.Neighbors, Verlet.Potentials, Verlet.Constraints, Verlet.Thermostats],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", "false") == "true"),
    checkdocs = :public,
        pages = [
        "Home" => "index.md",
        "Guide" => Any[
            "Constrained Dynamics" => "guide/constraints.md",
            "Forces & Potentials" => "guide/forces.md",
            "Bonded Potentials" => "guide/bonded.md",
            "Numerical Notes" => "guide/numerics.md",
        ],
        "Specs" => Any[
            "Overview" => "specs/Spec_Module_Index.md",
            "Top-level (Verlet)" => "specs/Spec_Verlet.md",
            "Core" => "specs/Spec_Core.md",
            "Neighbors" => "specs/Spec_Neighbors.md",
            "Potentials" => "specs/Spec_Potentials.md",
            "Constraints" => "specs/Spec_Constraints.md",
            "Thermostats" => "specs/Spec_Thermostats.md",
            "Electrostatics" => "specs/Spec_Electrostatics.md",
        ],
        "API" => "api.md",
    ],
)


deploydocs(; repo = "https://github.com/IlianPihlajamaa/Verlet.jl.git ", devbranch = "main")

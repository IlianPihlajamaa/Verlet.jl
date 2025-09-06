using Documenter
using Verlet

DocMeta.setdocmeta!(Verlet, :DocTestSetup, :(using Verlet); recursive=true)

makedocs(
    sitename = "Verlet.jl",
    modules = [Verlet],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", "false") == "true"),
    checkdocs = :public,
        pages = [
        "Home" => "index.md",
        "Guide" => Any[
            "Constrained Dynamics" => "guide/constraints.md",
            "Forces & Potentials" => "guide/forces.md",
            "Numerical Notes" => "guide/numerics.md",
        ],
        "API" => "api.md",
    ],
)


deploydocs(; repo = "https://github.com/IlianPihlajamaa/Verlet.jl.git ", devbranch = "main")
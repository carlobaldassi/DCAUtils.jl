using Documenter, DCAUtils

CIbuild = get(ENV, "CI", nothing) == "true"

makedocs(
    modules  = [DCAUtils],
    format   = Documenter.HTML(prettyurls = CIbuild),
    sitename = "DCAUtils.jl",
    pages    = Any[
        "Home" => "index.md",
       ]
    )

deploydocs(
    repo   = "github.com/carlobaldassi/DCAUtils.jl.git",
)

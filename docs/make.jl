using CalciumScoring
using Documenter

DocMeta.setdocmeta!(CalciumScoring, :DocTestSetup, :(using CalciumScoring); recursive=true)

makedocs(;
    modules=[CalciumScoring],
    authors="Dale <djblack@uci.edu> and contributors",
    repo="https://github.com/Dale-Black/CalciumScoring.jl/blob/{commit}{path}#{line}",
    sitename="CalciumScoring.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Dale-Black.github.io/CalciumScoring.jl",
        assets=String[],
    ),
    pages=["Home" => "index.md"],
)

deploydocs(; repo="github.com/Dale-Black/CalciumScoring.jl")

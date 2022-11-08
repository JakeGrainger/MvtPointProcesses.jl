using MvtPointProcesses
using Documenter

DocMeta.setdocmeta!(MvtPointProcesses, :DocTestSetup, :(using MvtPointProcesses); recursive=true)

makedocs(;
    modules=[MvtPointProcesses],
    authors="Jake Grainger",
    repo="https://github.com/JakeGrainger/MvtPointProcesses.jl/blob/{commit}{path}#{line}",
    sitename="MvtPointProcesses.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JakeGrainger.github.io/MvtPointProcesses.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JakeGrainger/MvtPointProcesses.jl",
    devbranch="main",
)

using GridMethod
using Documenter

DocMeta.setdocmeta!(GridMethod, :DocTestSetup, :(using GridMethod); recursive=true)

makedocs(;
    modules=[GridMethod],
    authors="Laura Brustenga i Moncusí <brust@math.ku.dk> and contributors",
    repo="https://github.com/LauraBMo/GridMethod.jl/blob/{commit}{path}#{line}",
    sitename="GridMethod.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://LauraBMo.github.io/GridMethod.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/LauraBMo/GridMethod.jl",
    devbranch="main",
)

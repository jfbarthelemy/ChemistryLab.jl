using CementChemistry
using Documenter

DocMeta.setdocmeta!(CementChemistry, :DocTestSetup, :(using CementChemistry); recursive=true)

makedocs(;
    modules=[CementChemistry],
    doctest = true,
    linkcheck = true,
    authors="Jean-François Barthélémy and Anthony Soive",
    repo = "https://jfbarthelemy.github.io/CementChemistry.jl/blob/{commit}{path}#{line}",
    sitename="CementChemistry.jl",
    format=Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical="https://jfbarthelemy.github.io/CementChemistry.jl",
        #edit_link="main",
        #assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Reference" => "reference.md",
        "Tutorial" => "tutorial.md",
    ],
)

deploydocs(;
    repo="github.com/jfbarthelemy/CementChemistry.jl",
    push_preview = false,
    devbranch="main",
)

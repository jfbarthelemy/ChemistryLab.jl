using CementChemistry
using Documenter

DocMeta.setdocmeta!(CementChemistry, :DocTestSetup, :(using CementChemistry); recursive=true)

makedocs(;
    modules=[CementChemistry],
    authors="Jean-François Barthélémy and Anthony Soive",
    sitename="CementChemistry.jl",
    format=Documenter.HTML(;
        canonical="https://jfbarthelemy.github.io/CementChemistry.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jfbarthelemy/CementChemistry.jl",
    devbranch="main",
)

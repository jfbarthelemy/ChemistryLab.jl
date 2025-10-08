using CementChemistry
using Documenter
using PrettyTables
using SymPy

DocMeta.setdocmeta!(CementChemistry, :DocTestSetup, :(using CementChemistry); recursive=true)

makedocs(;
    modules=[CementChemistry],
    # doctest = true,
    # linkcheck = true,
    authors="Jean-François Barthélémy and Anthony Soive",
    # repo = "https://github.com/jfbarthelemy/CementChemistry.jl/blob/{commit}{path}#{line}",
    sitename="CementChemistry.jl",
    format=Documenter.HTML(;
        canonical="https://jfbarthelemy.github.io/CementChemistry.jl",
        edit_link="main",
        assets = String[],
        ),
    pages=[
        "Home" => "index.md",
        "Tutorial" => [
            "man/formula_manipulation.md",
            "man/species.md",
            "man/databases.md",
            "man/stoichio_matrix.md",
            "man/equations.md",
            ],
        "Examples" => [
            "example/get_stochio_matrix.md",
            ],
        "API" => "reference.md",
    ],
)

deploydocs(;
    repo="github.com/jfbarthelemy/CementChemistry.jl",
    devbranch="main",
)

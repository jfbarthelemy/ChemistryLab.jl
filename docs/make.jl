using ChemistryLab
using Documenter
using PrettyTables

DocMeta.setdocmeta!(ChemistryLab, :DocTestSetup, :(using ChemistryLab); recursive=true)

makedocs(;
    modules=[ChemistryLab],
    # doctest = true,
    # linkcheck = true,
    authors="Jean-François Barthélémy and Anthony Soive",
    # repo = "https://github.com/jfbarthelemy/ChemistryLab.jl/blob/{commit}{path}#{line}",
    sitename="ChemistryLab.jl",
    format=Documenter.HTML(;
        canonical="https://jfbarthelemy.github.io/ChemistryLab.jl",
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
            "example/get_stoichio_matrix.md",
            "example/bogue_calculation.md",
            ],
        "API" => "reference.md",
    ],
)

deploydocs(;
    repo="github.com/jfbarthelemy/ChemistryLab.jl",
    devbranch="main",
)

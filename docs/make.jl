using CementChemistry
using Documenter

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
            "man/quickstart.md",
            "man/formula_manipulation.md",
            "man/databases.md",
            "man/stoichio_matrix.md",
            "tutorial.md",
            ],
        # "Tutorial" => "tutorial.md",
        "API" => "reference.md",
    ],
)

deploydocs(;
    repo="github.com/jfbarthelemy/CementChemistry.jl",
    devbranch="main",
)

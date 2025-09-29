push!(LOAD_PATH,"../src/")

using CementChemistry
using Documenter

DocMeta.setdocmeta!(CementChemistry, :DocTestSetup, :(using CementChemistry); recursive=true)

makedocs(;
    modules=[CementChemistry],
    # doctest = true,
    # linkcheck = true,
    # authors="Jean-François Barthélémy and Anthony Soive",
    # repo = "https://github.com/jfbarthelemy/CementChemistry.jl/blob/{commit}{path}#{line}",
    sitename="CementChemistry.jl",
    format=Documenter.HTML(;
        prettyurls = true,#get(ENV, "CI", "false") == "true",
        # canonical="https://jfbarthelemy.github.io/CementChemistry.jl",
        #edit_link="main",
        assets = ["assets/css/custom.css"],
        ),
    pages=[
        "Home" => "index.md",
        "Reference" => "reference.md",
        "Tutorial" => [
            # "man/quickstart.md",
            "tutorial.md",
            ],
        # "Tutorial" => "tutorial.md",
    ],
)

deploydocs(;
    repo="github.com/jfbarthelemy/CementChemistry.jl",
    push_preview = false,
    devbranch="main",
)

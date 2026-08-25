using ChemistryLab
using Optimization, OptimizationIpopt  # load extension OptimizationIpoptExt
using OptimaSolver                     # load extension OptimaSolverExt
using OrdinaryDiffEq                  # load extension KineticsOrdinaryDiffEqExt
using Documenter
using DocumenterCitations
using PrettyTables

include("pages.jl")

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style = :authoryear)

DocMeta.setdocmeta!(
    ChemistryLab,
    :DocTestSetup,
    :(using ChemistryLab, DynamicQuantities, OrderedCollections, Symbolics);
    recursive=true,
)

ENV["FORCE_COLOR"] = "true"
ENV["COLUMNS"] = "200"
ENV["LINES"] = "100"
ENV["GKSwstype"] = "100"   # headless GR backend — prevents Plots from hanging in doc builds

# Composite panels crop their outer labels unless the margins are explicit — a
# big enough canvas is not enough on its own. Setting the defaults once here
# covers every `@example` block, which run inside this process.
using Plots
Plots.gr()
Plots.default(;
    left_margin = 6Plots.mm,
    bottom_margin = 6Plots.mm,
    right_margin = 4Plots.mm,
    top_margin = 3Plots.mm,
)

makedocs(;
    clean=false,
    modules=[ChemistryLab],
    remotes=nothing,
    authors="Jean-François Barthélémy and Anthony Soive",
    sitename="ChemistryLab.jl",
    format=Documenter.HTML(;
        mathengine=Documenter.MathJax3(Dict(
            :loader => Dict("load" => ["[tex]/mhchem"]),
        )),
        canonical="https://MicroPoroChemoMechanics.github.io/ChemistryLab.jl",
        repolink="https://github.com/MicroPoroChemoMechanics/ChemistryLab.jl",
        edit_link="main",
        assets=["assets/favicon.ico", "assets/custom.css", "assets/mathjax-flash.js"],
        prettyurls=(get(ENV, "CI", nothing) == "true"),
        collapselevel=1,
        size_threshold_warn=200_000,
    ),
    pages=pages,
    plugins=[bib],
    warnonly=[:docs_block],
    draft=false,
)

deploydocs(;
    repo         = "github.com/MicroPoroChemoMechanics/ChemistryLab.jl.git",
    devbranch    = "main",
    push_preview = false,
)

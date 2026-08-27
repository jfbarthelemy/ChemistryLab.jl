using ChemistryLab
using Optimization, OptimizationIpopt  # load extension OptimizationIpoptExt
using OptimaSolver                     # load extension OptimaSolverExt
using OrdinaryDiffEq                  # load extension KineticsOrdinaryDiffEqExt
using Documenter
using DocumenterCitations
# VitePress renders the site from the Markdown that Documenter emits, and
# typesets every formula — chemical equations included — at build time into
# static SVG, so no MathJax bundle is ever fetched by the reader's browser.
using DocumenterVitepress
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

# ── Stopgap: CitationSiteNode in the Markdown writer ─────────────────────────
# DocumenterCitations 1.5 wraps every expanded citation in a `CitationSiteNode`,
# whose only purpose is to give the citation an HTML anchor so the bibliography
# can link back to it. Its own docstring calls it "transparent in any output
# format other than HTML", and both the LaTeX writer and MDFlatten implement it
# as "render my children".
#
# DocumenterVitepress 0.3.5 ships a DocumenterCitations extension, but it covers
# only `BibliographyNode`. With no method for `CitationSiteNode`, the writer
# falls through to its generic branch, which prints `Markdown.plain(element)` —
# so every one of this manual's citations came out as the literal text
# `DocumenterCitations.CitationSiteNode("kachanov1992-cite-1")`.
#
# The same one-line treatment as the other non-HTML writers. Remove this once
# DocumenterVitepress covers the node upstream.
function DocumenterVitepress.render(
        io::IO,
        mime::MIME"text/plain",
        node::Documenter.MarkdownAST.Node,
        ::DocumenterCitations.CitationSiteNode,
        page,
        doc;
        kwargs...,
    )
    return DocumenterVitepress.render(io, mime, node, node.children, page, doc; kwargs...)
end

# ── Stopgap: heading anchors that contain LaTeX ──────────────────────────────
# DocumenterVitepress builds each heading as `## <text> {#<slug>}`, where the
# slug is Documenter's anchor label passed through its own
# `sanitized_anchor_label` — whose comment says "vitepress doesn't like special
# markdown characters in the id slug", but which only strips `[ ] ( ) *`.
#
# A heading such as `## The COD tensor ``\boldsymbol{B}``` yields the slug
# `The-COD-tensor-\boldsymbol{B}`. VitePress's `{#...}` parser rejects the
# backslash and the braces, so it treats the whole suffix as *text*: the heading
# renders as "The COD tensor {#The-COD-tensor-\boldsymbol{B}}", the formula is
# dropped, and the same garbage lands in the "On this page" outline. Twenty-three
# headings across seven pages were affected.
#
# Stripping those characters from the slug is safe here: nothing links to those
# anchors (checked against every `](…#…)` destination in the generated Markdown),
# and this narrows to headings only, leaving docstring anchors — which legitimately
# carry braces, are emitted as raw `<a id=…>`, and *are* linked to — untouched.
#
# Remove once `sanitized_anchor_label` covers these characters upstream.
function DocumenterVitepress.render(
        io::IO,
        mime::MIME"text/plain",
        node::Documenter.MarkdownAST.Node,
        header::Documenter.AnchoredHeader,
        page,
        doc;
        kwargs...,
    )
    anchor = header.anchor
    label = DocumenterVitepress.sanitized_anchor_label(anchor)
    id = replace(replace(label, r"[\\{}]" => ""), " " => "-")
    heading = first(node.children)
    println(io)
    print(io, "#"^(heading.element.level), " ")
    heading_iob = IOBuffer()
    DocumenterVitepress.render(heading_iob, mime, node, heading.children, page, doc; kwargs...)
    print(io, rstrip(String(take!(heading_iob))))
    print(io, " {#$(id)}")
    if haskey(kwargs, :inventory)
        item = DocumenterVitepress.InventoryItem(
            name = anchor.id,
            domain = "std",
            role = "label",
            dispname = DocumenterVitepress._get_inventory_dispname(
                anchor.id, Documenter.MDFlatten.mdflatten(anchor.node)
            ),
            priority = -1,
            uri = DocumenterVitepress._get_inventory_uri(doc, page, id),
        )
        push!(kwargs[:inventory], item)
    end
    println(io)
    return nothing
end

# ── Stopgap: ordered lists start at 2, and swallow their first item ──────────
# DocumenterVitepress numbers ordered-list items with `bullet(i) = "$(i+1). "`,
# but `enumerate` is already 1-based: every ordered list in the manual came out
# numbered from 2. It also emits no blank line before the list.
#
# Together those two produce the damage seen on the References page. A list whose
# first marker is `2.` cannot interrupt a paragraph — CommonMark allows that only
# for a list starting at `1.` — so, with no blank line to separate them, the first
# entry was absorbed into the preceding prose as plain text and the list began at
# `3.`. Eighteen ordered lists across the manual were affected; not one of them
# started at 1.
#
# Remove once the numbering is fixed upstream.
function DocumenterVitepress.render(
        io::IO,
        mime::MIME"text/plain",
        node::Documenter.MarkdownAST.Node,
        list::Documenter.MarkdownAST.List,
        page,
        doc;
        kwargs...,
    )
    bullet(i) = list.type === :ordered ? "$(i). " : "- "
    println(io)
    iob = IOBuffer()
    for (i, item) in enumerate(node.children)
        DocumenterVitepress.render(
            iob, mime, item, item.children, page, doc; prenewline = false, kwargs...
        )
        eachline = split(String(take!(iob)), '\n')
        # Continuation lines must line up with the text, i.e. under the marker's
        # full width. Upstream hard-codes two spaces, which fits `- ` but not
        # `1. `: a display equation inside an ordered item fell out of the list,
        # splitting it in two and restarting the numbering.
        pad = " "^length(bullet(i))
        eachline[2:end] .= pad .* eachline[2:end]
        final_string = join(eachline, '\n')
        endswith(final_string, '\n') || (final_string *= "\n")
        print(io, bullet(i))
        print(io, final_string)
    end
    return nothing
end

makedocs(;
    # `clean = false` lets pages deleted from the source survive in `build/`
    # and go on being deployed. Nothing writes there before `makedocs`.
    modules=[ChemistryLab],
    remotes=nothing,
    authors="Jean-François Barthélémy and Anthony Soive",
    sitename="ChemistryLab.jl",
    # The favicon and the logo are picked up automatically from `docs/src/assets`,
    # and the sidebar is derived from `pages`, so neither needs declaring here.
    # mhchem now loads in `docs/src/.vitepress/mathjax-plugin.ts` instead of in a
    # `mathengine`, because the typesetting happens at build time.
    format=DocumenterVitepress.MarkdownVitepress(;
        repo = "https://github.com/MicroPoroChemoMechanics/ChemistryLab.jl",
        devbranch = "main",
        devurl = "dev",
        deploy_url = "https://MicroPoroChemoMechanics.github.io/ChemistryLab.jl",
        description = "Aqueous and cement chemistry in Julia: speciation, equilibrium and kinetics",
    ),
    pages=pages,
    plugins=[bib],
    warnonly=[:docs_block],
    draft=false,
)

# DocumenterVitepress writes a real directory per version rather than the
# symlinks Documenter used, so it needs its own `deploydocs`.
DocumenterVitepress.deploydocs(;
    repo         = "github.com/MicroPoroChemoMechanics/ChemistryLab.jl.git",
    target       = joinpath(@__DIR__, "build"),
    branch       = "gh-pages",
    devbranch    = "main",
    push_preview = false,
)

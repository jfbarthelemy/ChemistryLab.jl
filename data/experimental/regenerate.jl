# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)
#
# =============================================================================
#  regenerate.jl — rebuild the vendored calorimetry CSVs from their source
#
#  The two CSV files next to this script are subsets of a CC-BY-4.0 dataset
#  published on Zenodo. They are vendored so that the test suite and the
#  documentation build without network access, but nothing about them is
#  unverifiable: this script downloads the original deposit and rebuilds them
#  byte for byte.
#
#      Šmilauer V., Reiterman P. (2025)
#      Isothermal calorimetry database of 65 cements and approximations
#      Zenodo. https://doi.org/10.5281/zenodo.15212785
#      Licensed CC-BY-4.0.
#
#  Requires network access and the `unzip` executable. Neither is needed to
#  *use* the vendored files — only to regenerate them.
#
#  Usage:
#      julia data/experimental/regenerate.jl              # download and rebuild
#      julia data/experimental/regenerate.jl --check      # rebuild, compare, do
#                                                         # not overwrite
#      CALO_CACHE=/some/dir julia data/experimental/regenerate.jl
#                                                         # reuse a download
# =============================================================================

using Downloads
using Printf
using SHA

const ZENODO_RECORD = "15212785"
const ZENODO_DOI = "10.5281/zenodo.15212785"

const ZIP_MEASURED = "Isothermal calorimetry-cements.zip"
const ZIP_AFFINITY = "Fit-affinity.zip"

"""
Number of log-spaced sampling targets requested per record.

The originals hold ~12 000 rows on a *uniform* ~77 s grid, which weights the
ten-day plateau some hundred times more heavily than the acceleration peak. The
subset is drawn on a log-spaced grid instead, and by *nearest index* rather than
by interpolation, so every number written out is a measured value and not a
derived one. Consecutive targets that land on the same row collapse, so the
files come out somewhat shorter than this.
"""
const N_TARGET = 500

"""
The records to vendor: `(source stem, output stem, role)`.

`122` is the primary calibration target — its w/b of 0.50 and Blaine of
397 m²/kg are within a few percent of the CEM I 52.5 N formulation that
`scripts/ionic_hydration.jl` already runs after Lavergne et al. (2018), so the
one input the deposit does not report (the clinker phase composition) is the
least of an extrapolation it can be. `116` is the holdout: same nominal strength
class, different w/b *and* different Blaine, and a 26-day record.
"""
const RECORDS = [
    ("122-CEM I 52.5 R Cizkovice-397", "smilauer2025-122-cemI-52.5R-cizkovice", "calibration target"),
    ("116-CEM I 52.5 R Ladce-415", "smilauer2025-116-cemI-52.5R-ladce", "holdout"),
]

# ── download and extract ──────────────────────────────────────────────────────

"""
    zenodo_url(filename) -> String

Direct download URL of one file of the deposit.
"""
zenodo_url(filename) =
    "https://zenodo.org/records/$ZENODO_RECORD/files/$(replace(filename, " " => "%20"))?download=1"

"""
    fetch_deposit(cache) -> String

Download and unpack both archives of the deposit into `cache`, skipping whatever
is already there. Returns `cache`.
"""
function fetch_deposit(cache)
    mkpath(cache)
    for (zipname, subdir) in ((ZIP_MEASURED, "measured"), (ZIP_AFFINITY, "affinity"))
        dest = joinpath(cache, subdir)
        isdir(dest) && !isempty(readdir(dest)) && continue
        local_zip = joinpath(cache, zipname)
        if !isfile(local_zip)
            @info "downloading" zipname
            Downloads.download(zenodo_url(zipname), local_zip)
        end
        mkpath(dest)
        run(`unzip -o -q $local_zip -d $dest`)
    end
    return cache
end

# ── parsing the source format ─────────────────────────────────────────────────

"""
    parse_source(path) -> (header::Vector{String}, t, q, Q)

Split one file of the deposit into its quoted metadata lines and its three
numeric columns — hydration time (h), heat flow (W/g of binder), released heat
(J/g of binder).

The metadata block is neither fixed-length nor uniformly quoted across the
deposit, so the split is made on the first line that parses as three numbers
rather than on a line count.
"""
function parse_source(path)
    header = String[]
    t = Float64[]
    q = Float64[]
    Q = Float64[]
    for line in eachline(path)
        fields = split(strip(line))
        vals = length(fields) == 3 ? tryparse.(Float64, fields) : nothing
        if vals !== nothing && !any(isnothing, vals)
            push!(t, vals[1]); push!(q, vals[2]); push!(Q, vals[3])
        elseif isempty(t)
            push!(header, strip(line))
        end
    end
    isempty(t) && error("no numeric rows found in $path")
    return header, t, q, Q
end

"""
    header_value(header, key) -> Union{String, Nothing}

Value of a `"key: value"` metadata line, unquoted and stripped, or `nothing`.
"""
function header_value(header, key)
    for line in header
        clean = strip(line, ['"', ' '])
        startswith(clean, key * ":") || continue
        return strip(clean[(length(key) + 2):end], ['"', ' '])
    end
    return nothing
end

"""
    affinity_path(cache, stem) -> String

Path to the depositors' own affinity-model fit for this cement.
"""
affinity_path(cache, stem) =
    joinpath(cache, "affinity", "$(first(split(stem, '-')))-Fit-aff.csv")

"""
    affinity_line(cache, stem) -> Union{String, Nothing}

The authors' own four-parameter affinity-model fit for this cement, as the single
comment line their fit file carries. Kept for comparison; see the documentation
page.
"""
function affinity_line(cache, stem)
    path = affinity_path(cache, stem)
    isfile(path) || return nothing
    line = strip(first(eachline(path)))
    return startswith(line, "#") ? strip(lstrip(line, ['#', ' '])) : nothing
end

"""
    affinity_curve(cache, stem) -> Union{Tuple{Vector{Float64}, Vector{Float64}}, Nothing}

The fitted `(time [h], Q [J/g])` of the depositors' affinity model.

Their fit is carried along so that our residual can be read against somebody
else's on the very same record, instead of only against the database-wide figure
quoted in their paper. Reading their tabulated curve rather than reimplementing
their model is deliberate: a functional form recalled from memory is exactly the
kind of thing that is silently wrong.
"""
function affinity_curve(cache, stem)
    path = affinity_path(cache, stem)
    isfile(path) || return nothing
    t = Float64[]
    Q = Float64[]
    for line in eachline(path)
        startswith(strip(line), "#") && continue
        f = split(strip(line))
        length(f) >= 3 || continue
        v = tryparse.(Float64, f[1:3])
        any(isnothing, v) && continue
        push!(t, v[1]); push!(Q, v[3])
    end
    return isempty(t) ? nothing : (t, Q)
end

"""
    nearest(xs, x) -> Int

Index of the element of `xs` closest to `x`.
"""
nearest(xs, x) = argmin(abs.(xs .- x))

# ── subsetting ────────────────────────────────────────────────────────────────

"""
    log_subset(t, n) -> Vector{Int}

Indices of `t` closest to `n` log-spaced targets spanning `t`, deduplicated and
sorted, always keeping the first and last rows.

Nearest index, not interpolation: the point of the exercise is that every value
in the vendored file is one the calorimeter actually reported.
"""
function log_subset(t, n)
    targets = 10 .^ range(log10(first(t)), log10(last(t)); length = n)
    idx = [argmin(abs.(t .- τ)) for τ in targets]
    return sort!(unique!(vcat(1, idx, lastindex(t))))
end

# ── writing ───────────────────────────────────────────────────────────────────

"""
    render(source_stem, out_stem, role, cache) -> String

Full text of one vendored CSV: a `#`-commented provenance block, then a plain
comma-separated table.

The provenance block is not decoration. It carries the license the data comes
under, and it carries the three quantities a calibration must *not* fit because
the experiment measured them — the Blaine fineness, the w/b ratio and the
isothermal temperature.
"""
function render(source_stem, out_stem, role, cache)
    path = joinpath(cache, "measured", source_stem * ".csv")
    header, t, q, Q = parse_source(path)
    keep = log_subset(t, N_TARGET)
    aff = affinity_line(cache, source_stem)
    curve = affinity_curve(cache, source_stem)

    io = IOBuffer()
    println(io, "# Isothermal calorimetry of a CEM I paste — $role")
    println(io, "#")
    println(io, "# Source dataset: Šmilauer V., Reiterman P. (2025), Isothermal calorimetry")
    println(io, "#   database of 65 cements and approximations. Zenodo.")
    println(io, "#   doi: $ZENODO_DOI")
    println(io, "#   license: CC-BY-4.0  <https://creativecommons.org/licenses/by/4.0/>")
    println(io, "#   source file: $(source_stem).csv")
    println(io, "#")
    println(io, "# Companion article, to be cited (CC BY-NC-ND 3.0 — cite, do not reproduce):")
    println(io, "#   Šmilauer V., Edelmannová J., Reiterman P., Isothermal calorimetry database")
    println(io, "#   of 65 cements with analytical approximations, Ceramics-Silikáty 70(1),")
    println(io, "#   25-33 (2026). doi: 10.13168/cs.2025.0048")
    println(io, "#")
    for (key, label) in (
            "Cement name" => "cement", "Blaine" => "blaine", "w/b" => "wb",
            "Isothermal temperature" => "temperature",
        )
        v = header_value(header, key)
        v === nothing || println(io, "# $label: $v")
    end
    for line in header
        clean = strip(line, ['"', ' '])
        if startswith(clean, "Truncated from line") || startswith(clean, "Released heat up to")
            println(io, "# $clean")
        end
    end
    aff === nothing || println(io, "# reference fit by the depositors: $aff")
    println(io, "#")
    println(io, "# Subset of the source record: $(length(keep)) of $(length(t)) rows, chosen as")
    println(io, "#   the rows nearest to $N_TARGET log-spaced times. Values are verbatim; no")
    println(io, "#   interpolation, no smoothing. Regenerate with data/experimental/regenerate.jl")
    println(io, "#")
    if curve === nothing
        println(io, "# columns: hydration time (h), heat flow (W/g of binder), heat (J/g of binder)")
        println(io, "time_h,heat_flow_W_per_g,heat_J_per_g")
        for i in keep
            println(io, t[i], ",", q[i], ",", Q[i])
        end
    else
        println(io, "# The fourth column is NOT a measurement: it is the depositors' own")
        println(io, "#   affinity-model fit, read off their tabulated curve at the row nearest")
        println(io, "#   in time. It is carried so that a fit can be judged against somebody")
        println(io, "#   else's on the very same record.")
        println(io, "#")
        println(io, "# columns: hydration time (h), heat flow (W/g of binder), heat (J/g of")
        println(io, "#   binder), depositors' fitted heat (J/g of binder)")
        println(io, "time_h,heat_flow_W_per_g,heat_J_per_g,reference_fit_J_per_g")
        ta, Qa = curve
        for i in keep
            println(io, t[i], ",", q[i], ",", Q[i], ",", Qa[nearest(ta, t[i])])
        end
    end
    return String(take!(io))
end

# ── driver ────────────────────────────────────────────────────────────────────

function main(args = ARGS)
    check_only = "--check" in args
    cache = get(ENV, "CALO_CACHE", joinpath(tempdir(), "smilauer2025-zenodo-$ZENODO_RECORD"))
    fetch_deposit(cache)

    ok = true
    for (source_stem, out_stem, role) in RECORDS
        text = render(source_stem, out_stem, role, cache)
        out = joinpath(@__DIR__, out_stem * ".csv")
        digest = bytes2hex(sha256(text))
        if check_only
            same = isfile(out) && read(out, String) == text
            ok &= same
            @printf "%-46s %s  sha256 %s\n" out_stem (same ? "identical" : "DIFFERS  ") digest
        else
            write(out, text)
            @printf "%-46s %5d rows  sha256 %s\n" out_stem count(==('\n'), text) digest
        end
    end
    check_only && !ok && error("vendored files differ from what this script rebuilds")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

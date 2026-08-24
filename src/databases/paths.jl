# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

"""
    datapath(parts::AbstractString...) -> String

Absolute path to a file shipped in ChemistryLab's `data/` directory. Called
without argument, returns the directory itself.

This is the recommended way to name a bundled database, because it does not
depend on the working directory: a script written this way runs identically from
the package root, from an editor whose REPL started elsewhere, and inside a
documentation build.

# Examples

```julia
substances = build_species(datapath("cemdata18-thermofun.json"))
ss_phases  = build_solid_solutions(datapath("solid_solutions.toml"), dict)
readdir(datapath())                       # every bundled data file
datapath("experimental", "README.md")     # subdirectories work too
```

See also [`read_thermofun_database`](@ref), [`build_species`](@ref).
"""
function datapath(parts::AbstractString...)
    root = pkgdir(@__MODULE__)
    root === nothing && error(
        "cannot locate the ChemistryLab package directory, so the bundled " *
            "`data/` files cannot be resolved. Pass an explicit path instead.",
    )
    return joinpath(root, "data", parts...)
end

"""
    resolve_data_path(path::AbstractString) -> String

Resolve `path` to an existing file, trying in order:

 1. `path` as given — relative to the working directory, or absolute;
 2. `datapath(path)` — the same relative path under the bundled `data/`;
 3. `datapath(basename(path))` — a bundled data file of that name, whatever
    directory prefix the caller wrote;
 4. `joinpath(pkgdir(ChemistryLab), path)` — relative to the package root.

The working directory comes first, so a call that already resolves keeps
resolving to exactly the same file: the fallbacks can only turn a failure into a
success, never change an existing answer. And steps 2–3 can only succeed for a
name that *is* one of the bundled data files, so a mistyped path to a file of
one's own still fails loudly instead of silently loading something else.

Throws `ArgumentError` listing the bundled files when nothing matches.
"""
function resolve_data_path(path::AbstractString)
    isfile(path) && return String(path)

    for candidate in (datapath(path), datapath(basename(path)))
        isfile(candidate) && return candidate
    end

    root = pkgdir(@__MODULE__)
    if root !== nothing
        candidate = joinpath(root, path)
        isfile(candidate) && return candidate
    end

    bundled = try
        sort!(readdir(datapath()))
    catch
        String[]
    end
    throw(
        ArgumentError(
            "no such data file: \"$path\". It was looked for relative to the " *
                "working directory ($(pwd())), then among the data files shipped " *
                "with ChemistryLab ($(datapath())). Bundled files: " *
                "$(join(bundled, ", ")). Use `datapath(\"<name>\")` to name a " *
                "bundled file from anywhere.",
        ),
    )
end

"""
    display_data_path(path::AbstractString) -> String

Short, machine-independent label for `path`, meant for the banner a reader sees
rather than for opening a file.

A file inside the package is shown relative to the package root
(`data/cemdata18-thermofun.json`); anything else is shown unchanged. This
matters because these banners are captured verbatim into the documentation:
printing the absolute path would bake the build machine's directories
(`/home/runner/work/...`) into every page that loads a database.
"""
function display_data_path(path::AbstractString)
    root = pkgdir(@__MODULE__)
    if root !== nothing
        relative = relpath(abspath(path), root)
        startswith(relative, "..") || return relative
    end
    return String(path)
end

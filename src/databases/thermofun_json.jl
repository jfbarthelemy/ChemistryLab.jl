# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using DataFrames
using DynamicQuantities
using JSON
using OrderedCollections
using ProgressMeter
using Tables
using TOML

"""
    HKF_SI_CONVERSIONS

Hardcoded conversion factors from SUPCRT (cal, bar) to SI (J, Pa) for `eos_hkf_coeffs`.
The JSON unit metadata for a3 and a4 is incorrect (missing `/bar`), hence this explicit table.

| Symbol | SUPCRT unit             | SI unit             | Factor     |
|--------|-------------------------|---------------------|------------|
| a1     | cal/(mol·bar)           | J/(mol·Pa)          | 4.184e-5   |
| a2     | cal/mol                 | J/mol               | 4.184      |
| a3     | (cal·K)/(mol·bar)       | (J·K)/(mol·Pa)      | 4.184e-5   |
| a4     | cal·K/mol               | J·K/mol             | 4.184      |
| c1     | cal/(mol·K)             | J/(mol·K)           | 4.184      |
| c2     | cal·K/mol               | J·K/mol             | 4.184      |
| wref   | cal/mol                 | J/mol               | 4.184      |
"""
const HKF_SI_CONVERSIONS = OrderedDict{Symbol, Float64}(
    :a1 => 4.184e-5,
    :a2 => 4.184,
    :a3 => 4.184e-5,
    :a4 => 4.184,
    :c1 => 4.184,
    :c2 => 4.184,
    :wref => 4.184,
)

"""
    read_thermofun_database(filename::AbstractString) -> (DataFrame, DataFrame, DataFrame)

Read a ThermoFun database from a JSON file.

# Arguments

  - `filename`: path to the JSON database file.

# Returns

  - `df_elements`: DataFrame of chemical elements.
  - `df_substances`: DataFrame of chemical substances (species).
  - `df_reactions`: DataFrame of chemical reactions.
"""
function read_thermofun_database(filename)
    path = resolve_data_path(filename)
    print_title(
        # `display_data_path`, not `path`: the banner is captured into the
        # documentation, so it must not carry a build machine's directories.
        "Loading database: $(display_data_path(path))";
        crayon = Crayon(; foreground = :green),
        style = :box,
        indent = "",
    )
    data = JSON.parsefile(path)
    df_substances = DataFrame(Tables.dictrowtable(data["substances"]))
    df_reactions = DataFrame(Tables.dictrowtable(data["reactions"]))
    df_elements = DataFrame(Tables.dictrowtable(data["elements"]))
    return df_elements, df_substances, df_reactions
end

"""
    extract_unit(v, default_unit=u"1") -> AbstractQuantity

Try to parse the string `v` as a unit via `uparse`.
Returns `default_unit` if parsing fails.
"""
function extract_unit(v, default_unit = u"1")
    return try
        uparse(v)
    catch
        default_unit
    end
end

"""
    extract_value(row, field; verbose=false, default_unit=u"1", with_units=true) -> Union{AbstractQuantity, Number, Missing}

Extract a scalar value (optionally with units) from a nested ThermoFun DataFrame `row`.
Returns `missing` when the field is absent, missing, or cannot be parsed.
"""
function extract_value(
        row, field::Symbol; verbose = false, default_unit = u"1", with_units = true
    )
    if haskey(row, field) && !ismissing(row[field]) && haskey(row[field], :values)
        try
            val = only(row[field].values)
            if with_units
                if iszero(val)
                    val *= default_unit
                elseif haskey(row[field], :units)
                    vunit = only(get(row[field], :units, [""]))
                    val *= extract_unit(vunit, default_unit)
                else
                    val *= default_unit
                end
            end
            if verbose
                println("$(row.symbol) => $field=$val")
            end
            return val
        catch
            return missing
        end
    else
        return missing
    end
end

correct_volume_unit(v::AbstractQuantity) = uamount(v) != -1 ? v / 1u"mol" : v

correct_volume_unit(v) = v

"""
    complete_species_with_thermo_model!(species, row; verbose=false)

Populate thermodynamic reference values and build thermodynamic functions on `species`
from a ThermoFun substance DataFrame `row`. Mutates `species.properties` in place.
"""
function complete_species_with_thermo_model!(species, row; verbose = false)
    Tref = row.Tst * u"K"
    Pref = row.Pst * u"Pa"
    species.Tref = Tref
    species.Pref = Pref
    values0 = [
        :Cp⁰ => extract_value(
            row, :sm_heat_capacity_p; verbose = verbose, default_unit = u"J/K/mol"
        ),
        :ΔₐH⁰ => extract_value(row, :sm_enthalpy; verbose = verbose, default_unit = u"J/mol"),
        :S⁰ =>
            extract_value(row, :sm_entropy_abs; verbose = verbose, default_unit = u"J/K/mol"),
        :ΔₐG⁰ =>
            extract_value(row, :sm_gibbs_energy; verbose = verbose, default_unit = u"J/mol"),
        :V⁰ => correct_volume_unit(extract_value(row, :sm_volume; verbose = verbose, default_unit = u"J/bar")),
    ]
    species[:thermo_params] = [values0; :T => Tref; :P => Pref]
    TPMethods = row.TPMethods
    if !ismissing(TPMethods)
        for method in TPMethods
            method_type = only(values(method.method))
            if method_type == "cp_ft_equation" && haskey(method, :m_heat_capacity_ft_coeffs)
                species[:thermo_method] = "cp_ft_equation"
                coeffs = method.m_heat_capacity_ft_coeffs
                vals = coeffs.values
                units = extract_unit.(coeffs.units)
                params = [
                    Symbol("a", subscriptnumber(i - 1)) => float(vals[i] * units[i]) for
                        i in 1:min(length(vals), length(units))
                ]
                species[:thermo_params] = [params; species[:thermo_params]]

            elseif method_type == "solute_hkf88_reaktoro" && haskey(method, :eos_hkf_coeffs)
                species[:thermo_method] = "solute_hkf88_reaktoro"
                coeffs = method.eos_hkf_coeffs
                vals = float.(coeffs.values)
                names = [:a1, :a2, :a3, :a4, :c1, :c2, :wref]
                hkf_params = [
                    names[i] => vals[i] * HKF_SI_CONVERSIONS[names[i]] for
                        i in 1:min(length(vals), length(names))
                ]
                z = float(get(row, :formula_charge, 0))
                push!(hkf_params, :z => z)
                species[:thermo_params] = [hkf_params; species[:thermo_params]]

            elseif method_type == "mv_constant"
                species[:V_method] = "mv_constant"
            end
        end
    end
    return species
end

"""
    build_species(df_substances::AbstractDataFrame, list_symbols=nothing; verbose=false) -> Vector{Species}

Build Species objects from a substance DataFrame.

# Arguments

  - `df_substances`: DataFrame containing substance data.
  - `list_symbols`: optional list of symbols to filter (default: nothing, process all).
  - `verbose`: if true, print details during processing (default: false).

# Returns

  - Vector of `Species`.
"""
function build_species(
        df_substances::AbstractDataFrame, list_symbols = nothing; verbose = false
    )
    local_df_substances = if isnothing(list_symbols)
        df_substances
    else
        @view df_substances[df_substances.symbol .∈ Ref(list_symbols), :]
    end
    keylist = String[]
    species_list = Species[]
    print_title(
        "Building species"; crayon = Crayon(; foreground = :blue), style = :box, indent = ""
    )
    @showprogress for row in eachrow(local_df_substances)
        if verbose
            println(row[:symbol])
        end
        species = Species(
            row.formula;
            name = row.name,
            symbol = row.symbol,
            aggregate_state = try
                eval(Meta.parse(only(values(row.aggregate_state))))
            catch
                AS_UNDEF
            end,
            class = try
                eval(Meta.parse(only(values(row.class_))))
            catch
                SC_UNDEF
            end,
        )
        complete_species_with_thermo_model!(species, row; verbose = verbose)
        key = row.symbol
        if key in keylist
            @warn("Symbol $key is used for multiple species")
        end
        push!(keylist, key)
        push!(species_list, species)
    end
    return species_list
end

function build_species(filename, list_symbols = nothing; verbose = false)
    _, df_substances, _ = read_thermofun_database(filename)
    return build_species(df_substances, list_symbols; verbose = verbose)
end

"""
    complete_reaction_with_thermo_model!(reaction, row; verbose=false)

Populate thermodynamic reference values and build thermodynamic functions on `reaction`
from a ThermoFun reaction DataFrame `row`. Mutates `reaction.properties` in place.
"""
function complete_reaction_with_thermo_model!(reaction, row; verbose = false)
    Tref = row.Tst * u"K"
    Pref = row.Pst * u"Pa"
    reaction.Tref = Tref
    reaction.Pref = Pref
    values0 = [
        :ΔᵣCp⁰ => extract_value(
            row, :drsm_heat_capacity_p; verbose = verbose, default_unit = u"J/K/mol"
        ),
        :ΔᵣH⁰ => extract_value(row, :drsm_enthalpy; verbose = verbose, default_unit = u"J/mol"),
        :ΔᵣS⁰ =>
            extract_value(row, :drsm_entropy_abs; verbose = verbose, default_unit = u"J/K/mol"),
        :ΔᵣG⁰ =>
            extract_value(row, :drsm_gibbs_energy; verbose = verbose, default_unit = u"J/mol"),
        :ΔᵣV⁰ => correct_volume_unit(extract_value(row, :drsm_volume; verbose = verbose, default_unit = u"J/bar")),
        :logKr => extract_value(row, :logKr; verbose = verbose, default_unit = u"1"),
    ]
    reaction[:thermo_params] = [values0; :T => Tref; :P => Pref]
    TPMethods = row.TPMethods
    if !ismissing(TPMethods)
        for method in TPMethods
            method_type = only(values(method.method))
            if method_type == "logk_fpt_function" && haskey(method, :logk_ft_coeffs)
                reaction[:logk_method] = "logk_fpt_function"
                coeffs = method.logk_ft_coeffs
                vals = coeffs.values
                units = dimension.([1, u"1/K", u"K", 1, u"K^2", u"1/K^2", u"1/√K"])
                params = [
                    Symbol("A", subscriptnumber(i - 1)) =>
                        float(Quantity(vals[i], units[i])) for
                        i in 1:min(length(vals), length(units))
                ]
                reaction[:thermo_params] = [params; reaction[:thermo_params]]
            elseif method_type == "dr_volume_constant"
                reaction[:V_method] = "dr_volume_constant"
            end
        end
    end
    return reaction
end

"""
    build_reactions(df_reactions::AbstractDataFrame, dict_species=Dict(), list_symbols=nothing; verbose=false) -> Vector{Reaction}

Build Reaction objects from a reaction DataFrame.

# Arguments

  - `df_reactions`: DataFrame containing reaction data.
  - `species_list`: vector of existing `Species` objects to use in reactions.
  - `list_symbols`: optional list of reaction symbols to filter (default: nothing, process all).
  - `verbose`: if true, print details during processing (default: false).

# Returns

  - Vector of `Reaction` objects.
"""
function build_reactions(
        df_reactions::AbstractDataFrame,
        species_list = [],
        list_symbols = nothing;
        verbose = false,
    )
    local_df_reactions = if isnothing(list_symbols)
        df_reactions
    else
        @view df_reactions[df_reactions.symbol .∈ Ref(list_symbols), :]
    end
    dict_species = Dict(symbol(s) => s for s in species_list)
    keylist = String[]
    reactions_list = Reaction[]
    print_title(
        "Building reactions"; crayon = Crayon(; foreground = :red), style = :box, indent = ""
    )
    function choose_species(k, rowsymbol, dict_species)
        if haskey(dict_species, k)
            return dict_species[k]
        elseif haskey(dict_species, rowsymbol)
            return dict_species[rowsymbol]
        else
            rowsymboldot = replace(rowsymbol, "_" => ".")
            if haskey(dict_species, rowsymboldot)
                return dict_species[rowsymboldot]
            else
                return find_species(k, collect(values(dict_species)))
            end
        end
    end
    @showprogress for row in eachrow(local_df_reactions)
        if verbose
            println(row[:symbol])
        end
        reaction = Reaction(
            OrderedDict(
                choose_species(last(k), row.symbol, dict_species) => last(v) for
                    (k, v) in row.reactants if last(k) != "e-"
            );
            symbol = row.symbol,
        )
        complete_reaction_with_thermo_model!(reaction, row; verbose = verbose)
        key = row.symbol
        if key in keylist
            @warn("Symbol $key is used for multiple reactions")
        end
        push!(keylist, key)
        push!(reactions_list, reaction)
    end
    return reactions_list
end

function build_reactions(filename, species_list = [], list_symbols = nothing; verbose = false)
    _, _, df_reactions = read_thermofun_database(filename)
    return build_reactions(df_reactions, species_list, list_symbols; verbose = verbose)
end

"""
    get_compatible_species(df_substances::AbstractDataFrame, species_list; aggregate_states=[AS_AQUEOUS], exclude_species=[], union=false) -> DataFrame

Find species in the database compatible with a given list of species (sharing atoms).

# Arguments

  - `df_substances`: substance DataFrame.
  - `species_list`: list of target species symbols.
  - `aggregate_states`: filter for specific aggregate states (default: `[AS_AQUEOUS]`).
  - `exclude_species`: list of species symbols to exclude.
  - `union`: if true, includes the original `species_list` in the result (default: false).

# Returns

  - DataFrame of compatible substances.
"""
function get_compatible_species(
        df_substances::AbstractDataFrame,
        species_list;
        aggregate_states = [AS_AQUEOUS],
        exclude_species = [],
        union = false,
    )
    df_given_species = @view df_substances[df_substances.symbol .∈ Ref(species_list), :]
    involved_atoms = union_atoms(parse_formula.(df_given_species.formula))
    mask1 = last.(only.(df_substances.aggregate_state)) .∈ Ref(string.(aggregate_states))
    mask2 = issubset.(keys.(parse_formula.(df_substances.formula)), Ref(involved_atoms))
    mask3 = .!(df_substances.symbol .∈ Ref(exclude_species))
    df_compat = @view df_substances[mask1 .&& mask2 .&& mask3, :]
    if union
        return unique(vcat(df_given_species, df_compat))
    else
        return df_compat
    end
end

"""
    build_solid_solutions(toml_file, dict_species; skip_missing=true) -> Vector{SolidSolutionPhase}

Load solid solution phase definitions from a TOML file and assemble
[`SolidSolutionPhase`](@ref) objects from an existing species dictionary.

Each end-member species is automatically requalified to `SC_SSENDMEMBER` via
[`with_class`](@ref), regardless of the class stored in the database.

# Arguments

  - `toml_file`: path to a TOML file with `[[solid_solution]]` entries (see
    `data/solid_solutions.toml` for the format).
  - `dict_species`: `Dict{String, <:AbstractSpecies}` mapping symbol → species
    (typically built from `Dict(symbol(s) => s for s in build_species(...))`).
  - `skip_missing`: if `true` (default), silently skip phases whose end-members
    are not all present in `dict_species`; if `false`, throw an error.

# TOML format

```toml
[[solid_solution]]
name        = "CSHQ"
end_members = ["CSHQ-TobD", "CSHQ-TobH", "CSHQ-JenH", "CSHQ-JenD"]
model       = "ideal"          # or "redlich_kister"
# For redlich_kister only:
a0          = 3000.0           # J/mol
a1          = 500.0            # J/mol
a2          = 0.0              # J/mol
```

# Example

```julia
substances = build_species(datapath("cemdata18-thermofun.json"))
dict       = Dict(symbol(s) => s for s in substances)
ss_phases  = build_solid_solutions(datapath("solid_solutions.toml"), dict)
cs = ChemicalSystem(species, CEMDATA_PRIMARIES; solid_solutions = ss_phases)
```
"""
function build_solid_solutions(
        toml_file::AbstractString,
        dict_species::AbstractDict;
        skip_missing::Bool = true,
    )
    data = TOML.parsefile(resolve_data_path(toml_file))
    entries = get(data, "solid_solution", [])
    phases = SolidSolutionPhase[]
    for entry in entries
        ss_name = entry["name"]
        em_symbols = entry["end_members"]

        # Check all end-members are available
        missing_syms = filter(sym -> !haskey(dict_species, sym), em_symbols)
        if !isempty(missing_syms)
            missing_str = join(missing_syms, ", ")
            if skip_missing
                @warn "build_solid_solutions: skipping \"$ss_name\" — " *
                    "end-members not found in dict_species: $missing_str"
                continue
            else
                error(
                    "build_solid_solutions: end-members not found for \"$ss_name\": " *
                        missing_str,
                )
            end
        end

        em_species = [dict_species[sym] for sym in em_symbols]

        # Build mixing model
        model_str = get(entry, "model", "ideal")
        mixing_model = if model_str == "redlich_kister"
            RedlichKisterModel(;
                a0 = get(entry, "a0", 0.0),
                a1 = get(entry, "a1", 0.0),
                a2 = get(entry, "a2", 0.0),
            )
        elseif model_str == "ideal"
            IdealSolidSolutionModel()
        else
            @warn "build_solid_solutions: unknown model \"$model_str\" for " *
                "\"$ss_name\", defaulting to ideal"
            IdealSolidSolutionModel()
        end

        push!(phases, SolidSolutionPhase(ss_name, em_species; model = mixing_model))
    end
    return phases
end

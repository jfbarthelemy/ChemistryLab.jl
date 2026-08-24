# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using JSON
using TOML

"""
    parse_reaction_stoich_cemdata(reaction_line::AbstractString) -> (Vector, String, String)

Parse a reaction line from CEMDATA format and extract stoichiometric information.

# Arguments

  - `reaction_line`: reaction string in CEMDATA format, optionally with a comment after '#'.

# Returns

  - `reactants`: vector of dictionaries with "symbol" and "coefficient" keys.
  - `modified_equation`: equation string with added "@" markers for aqueous species.
  - `comment`: extracted comment string (empty if none).

The function automatically adds "@" suffixes to aqueous species without explicit charges
(except for the first reactant).
"""
function parse_reaction_stoich_cemdata(reaction_line::AbstractString)
    # Split the equation and the comment
    equation_parts = split(reaction_line, '#')
    equation = strip(equation_parts[1])
    comment = length(equation_parts) > 1 ? strip(equation_parts[2]) : ""

    parts = split(equation, "=")
    if length(parts) != 2
        @warn "Equation format is incorrect: $equation"
        return [], equation, comment
    end

    modified_equation_parts = String[]
    reactants = []

    for (side_index, side) in enumerate(parts)
        tokens = split(side)
        modified_tokens = []
        for (token_index, token) in enumerate(tokens)
            if token == "+"
                push!(modified_tokens, token)
                continue
            end
            m = match(r"^([+-]?\d*\.?\d+)?(.+?)([\+\-]\d*\.?\d*)?$", token)
            if m !== nothing
                coeff_str = m.captures[1]
                base_symbol = m.captures[2]
                charge = m.captures[3]
                if isnothing(charge)
                    charge = ""
                end
                sp = base_symbol * charge
                coeff = if isnothing(coeff_str) || isempty(coeff_str)
                    1.0
                else
                    parse(Float64, coeff_str)
                end

                # Add "@" to aqueous species without charge
                if isempty(charge) && !(side_index == 1 && token_index == 1)
                    sp *= "@"
                    token = coeff_str === nothing ? sp : coeff_str * sp
                end

                push!(modified_tokens, token)

                # Add to reactants/products
                push!(
                    reactants,
                    Dict("symbol" => sp, "coefficient" => side_index == 1 ? -coeff : coeff),
                )
            else
                push!(modified_tokens, token)
            end
        end
        push!(modified_equation_parts, join(modified_tokens, " "))
    end
    modified_equation = join(modified_equation_parts, " = ")

    return reactants, modified_equation, comment
end

"""
    parse_float_array(line::AbstractString) -> Vector{Float64}

Parse a line containing space-separated floats, skipping the first token and any comments.

# Arguments

  - `line`: input string with format "keyword value1 value2 ...".

# Returns

  - Vector of successfully parsed Float64 values.

# Examples

```julia
julia> parse_float_array("-analytical_expression 1.5 2.3 4.7")
3-element Vector{Float64}:
 1.5
 2.3
 4.7

julia> parse_float_array("-log_K 5.2 # comment")
1-element Vector{Float64}:
 5.2
```
"""
function parse_float_array(line)
    parts = split(line)
    if length(parts) < 2
        return Float64[]
    end
    float_parts = Float64[]
    for part in parts[2:end]
        if !startswith(part, "#")
            try
                push!(float_parts, parse(Float64, part))
            catch e
                # @warn "Could not parse '$part' as Float64, skipping."
            end
        end
    end
    return float_parts
end

"""
    parse_phases(dat_content::AbstractString) -> Dict{String,Any}

Extract phase information from PHREEQC .dat file content.

# Arguments

  - `dat_content`: full text content of a PHREEQC .dat file.

# Returns

  - Dictionary mapping phase names to their properties (equation, log_K, analytical_expression, V⁰).

Parses the PHASES section and extracts reaction equations, equilibrium constants, analytical
expressions, and molar volumes for each phase.
"""
function parse_phases(dat_content)
    phases = Dict{String, Any}()
    in_phases = false
    current_phase = nothing

    for line in eachline(IOBuffer(dat_content))
        line = strip(line)
        if startswith(line, "PHASES")
            in_phases = true
            continue
        elseif in_phases && !isempty(line) && !startswith(line, "#")
            if !occursin("=", line) && !startswith(line, "-") && !isempty(line)
                parts = split(line)
                if length(parts) >= 1 && !startswith(parts[1], "-")
                    phase_name = parts[1]
                    current_phase = Dict{String, Any}("symbol" => phase_name)
                    phases[phase_name] = current_phase
                end
            elseif occursin("=", line) && current_phase !== nothing
                reactants, equation, comment = parse_reaction_stoich_cemdata(line)
                current_phase["equation"] = equation
                current_phase["reactants"] = reactants
                if !isempty(comment)
                    current_phase["comment"] = comment
                end
            elseif startswith(line, "-log_K") && current_phase !== nothing
                log_k_parts = split(line)
                if length(log_k_parts) >= 2
                    try
                        current_phase["logKr"] = Dict(
                            "values" => [parse(Float64, log_k_parts[2])], "errors" => [2]
                        )
                    catch e
                        @warn "Could not parse log_K value for phase $(current_phase["symbol"]), skipping."
                    end
                end
            elseif startswith(line, "-analytical_expression") && current_phase !== nothing
                analytical_expression = parse_float_array(line)
                # coef of log10 in .dat becomes a coef of log in .json
                if length(analytical_expression) > 3
                    analytical_expression[4] /= log(10)
                end
                current_phase["analytical_expression"] = analytical_expression
            elseif startswith(line, "-V⁰") && current_phase !== nothing
                V⁰_parts = split(line)
                if length(V⁰_parts) >= 2
                    try
                        current_phase["drsm_volume"] = parse(Float64, V⁰_parts[2])
                    catch e
                        @warn "Could not parse V⁰ value for phase $(current_phase["symbol"]), skipping."
                    end
                end
            end
        end
    end

    return phases
end

"""
    extract_primary_species(file_path::AbstractString) -> DataFrame

Extract primary aqueous species from a PHREEQC database file.

# Arguments

  - `file_path`: path to PHREEQC .dat file.

# Returns

  - DataFrame with columns: species, symbol, formula, aggregate_state, atoms, charge, gamma.

Parses the SOLUTION_SPECIES section to extract master species and their properties.
The "Zz" charge placeholder is handled specially. Gamma coefficients for activity
models are extracted from "-gamma" lines.
"""
function extract_primary_species(file_path)
    lines = readlines(resolve_data_path(file_path))

    start_idx = 0
    end_idx = 0
    in_primary_section = false

    for (i, line) in enumerate(lines)
        stripped = strip(line)

        if startswith(stripped, "SOLUTION_SPECIES")
            start_idx = i + 1
            in_primary_section = true

            while start_idx <= length(lines)
                next_line = strip(lines[start_idx])
                if startswith(next_line, "# PMATCH MASTER SPECIES") ||
                        occursin("=", next_line)
                    break
                end
                start_idx += 1
            end

            if startswith(strip(lines[start_idx]), "# PMATCH MASTER SPECIES")
                start_idx += 1
            end
        end

        if in_primary_section && startswith(stripped, "# PMATCH SECONDARY MASTER SPECIES")
            end_idx = i - 1
            break
        end
    end

    if in_primary_section && end_idx == 0
        end_idx = length(lines)
    end

    species_data = []

    for i in start_idx:end_idx
        line = strip(lines[i])

        if occursin("=", line)
            parts = split(line, "=")
            current_species = strip(parts[1])

            if current_species == "e-"
                symbol = "Zz"
            elseif occursin(r"[\+\-]\d*$", current_species)
                symbol = current_species
            else
                symbol = current_species * "@"
            end

            push!(species_data, (species = current_species, symbol = symbol, gamma = Float64[]))
        end

        if startswith(line, "-gamma") && !isempty(species_data)
            parts = split(line)
            gamma_values = Float64[]
            for val in parts[2:end]
                num = tryparse(Float64, val)
                if num !== nothing
                    push!(gamma_values, num)
                end
            end

            if !isempty(gamma_values)
                last_entry = species_data[end]
                species_data[end] = (
                    species = last_entry.species, symbol = last_entry.symbol, gamma = gamma_values,
                )
            end
        end
    end

    df = DataFrame(species_data)
    df.symbol = String.(df.symbol)
    df.formula .= df.symbol
    df.aggregate_state .= "AS_AQUEOUS"
    df.atoms .= parse_formula.(df.symbol)
    df.charge .= extract_charge.(df.symbol)
    df[df.symbol .== "Zz", :species] .= "Zz"
    df[df.symbol .== "Zz", :formula] .= "Zz"
    df[df.symbol .== "Zz", :charge] .= 1
    return df[sortperm(df.symbol .== "Zz"), :]
end

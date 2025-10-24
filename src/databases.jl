"""
    get_aqueous_species(json_data)

Extract the set of aqueous species from a ThermoFun JSON data structure.

# Arguments
- `json_data`: Parsed JSON object containing a "substances" field.

# Returns
- `Set{String}`: Set of aqueous species symbols.
"""
function get_aqueous_species(json_data)
    aqueous_species = Set{String}()
    for substance in json_data["substances"]
        aggregate_state = get(substance, "aggregate_state", nothing)
        if aggregate_state isa Dict && haskey(aggregate_state, "4") && aggregate_state["4"] == "AS_AQUEOUS"
            str = substance["symbol"]
            if str[end] == '@'
                str = first(str, length(str)-1)
            end
            push!(aqueous_species, str)
        end
    end
    return aqueous_species
end

"""
    parse_reaction_stoich_cemdata(reaction_line::AbstractString, aqueous_species, gaseous=false)

Parse a Cemdata/Phreeqc reaction line, adding "@" for aqueous species as needed.

# Arguments
- `reaction_line::AbstractString`: The reaction line from a Cemdata .dat file.
- `aqueous_species`: Set of aqueous species symbols.
- `gaseous`: Boolean indicating if the phase is gaseous.

# Returns
- `reactants`: Array of Dicts with "symbol" and "coefficient" for each reactant/product.
- `modified_equation`: The formatted equation string.
- `comment`: Any comment found on the line.
"""
function parse_reaction_stoich_cemdata(reaction_line::AbstractString, aqueous_species, gaseous=false)
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
            m = match(r"^([+-]?\d*\.?\d+)?([A-Za-z][A-Za-z0-9\(\)]*)([\+\-]?\d*)$", token)
            if m !== nothing
                coeff_str = m.captures[1]
                base_symbol = m.captures[2]
                charge = m.captures[3]
                sp = base_symbol * charge
                coeff = coeff_str === nothing || isempty(coeff_str) ? 1.0 : parse(Float64, coeff_str)

                # Add "@" to aqueous species without charge
                if base_symbol in aqueous_species && isempty(charge) && !(gaseous && side_index == 1 && token_index == 1)
                    sp *= "@"
                    token = coeff_str === nothing ? sp : coeff_str * sp
                end

                push!(modified_tokens, token)

                # Add to reactants/products
                push!(reactants, Dict("symbol" => sp, "coefficient" => side_index == 1 ? -coeff : coeff))
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
    parse_float_array(line)

Parse a line containing a float array (e.g., -analytical_expression) and return the array of Float64 values.

# Arguments
- `line`: The line to parse.

# Returns
- `Vector{Float64}`: Array of parsed float values.
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
                @warn "Could not parse '$part' as Float64, skipping."
            end
        end
    end
    return float_parts
end

"""
    parse_phases(dat_content, aqueous_species)

Parse the PHASES section of a Cemdata .dat file, extracting phase info and reactions.

# Arguments
- `dat_content`: The content of the Cemdata .dat file as a string.
- `aqueous_species`: Set of aqueous species symbols.

# Returns
- `Dict{String, Any}`: Dictionary of phase names to their data.
"""
function parse_phases(dat_content, aqueous_species)
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
                reactants, equation, comment = parse_reaction_stoich_cemdata(line, aqueous_species, occursin("(g)", get(current_phase, "symbol", "")))
                current_phase["equation"] = equation
                current_phase["reactants"] = reactants
                if !isempty(comment)
                    current_phase["comment"] = comment
                end
            elseif startswith(line, "-log_K") && current_phase !== nothing
                log_k_parts = split(line)
                if length(log_k_parts) >= 2
                    try
                        current_phase["logKr"] = Dict("values" => [parse(Float64, log_k_parts[2])], "errors" => [2])
                    catch e
                        @warn "Could not parse log_K value for phase $(current_phase["symbol"]), skipping."
                    end
                end
            elseif startswith(line, "-analytical_expression") && current_phase !== nothing
                current_phase["analytical_expression"] = parse_float_array(line)
            elseif startswith(line, "-Vm") && current_phase !== nothing
                vm_parts = split(line)
                if length(vm_parts) >= 2
                    try
                        current_phase["drsm_volume"] = parse(Float64, vm_parts[2])
                    catch e
                        @warn "Could not parse Vm value for phase $(current_phase["symbol"]), skipping."
                    end
                end
            end
        end
    end

    return phases
end

"""
    merge_reactions(json_data, new_reactions)

Merge new reactions into a ThermoFun JSON structure, skipping duplicates.

# Arguments
- `json_data`: The parsed JSON object containing a "reactions" field.
- `new_reactions`: Dictionary of new reactions to add.

# Returns
- The updated JSON object with merged reactions.
"""
function merge_reactions(json_data, new_reactions)
    existing_symbols = Set{String}()
    for reaction in json_data["reactions"]
        push!(existing_symbols, reaction["symbol"])
    end

    new_reactions_list = []
    for (name, phase) in new_reactions
        if !(name in existing_symbols)
            if haskey(phase, "logKr") && haskey(phase, "analytical_expression") && haskey(phase, "equation")
                reaction_dict = Dict{String, Any}()

                reaction_dict["symbol"] = phase["symbol"]
                reaction_dict["equation"] = phase["equation"]
                if haskey(phase, "comment")
                    reaction_dict["comment"] = phase["comment"]
                end
                reaction_dict["reactants"] = phase["reactants"]

                reaction_dict["limitsTP"] = Dict{String, Any}(
                    "range" => false,
                    "lowerP" => 0.1,
                    "lowerT" => 273.15,
                    "upperP" => 1000000,
                    "upperT" => 298.15
                )

                reaction_dict["Tst"] = 298.15
                reaction_dict["Pst"] = 100000

                reaction_dict["TPMethods"] = [
                    Dict{String, Any}(
                        "method" => Dict("0" => "logk_fpt_function"),
                        "limitsTP" => Dict{String, Any}(
                            "lowerP" => 0,
                            "lowerT" => 273.15,
                            "upperP" => 0,
                            "upperT" => 273.15
                        ),
                        "logk_ft_coeffs" => Dict{String, Any}(
                            "values" => vcat(phase["analytical_expression"], zeros(max(0, 12 - length(phase["analytical_expression"]))))
                        )
                    ),
                    Dict{String, Any}(
                        "method" => Dict("7" => "logk_3_term_extrap")
                    ),
                    Dict{String, Any}(
                        "method" => Dict("13" => "dr_volume_constant")
                    )
                ]

                reaction_dict["logKr"] = phase["logKr"]

                R = 8.31446261815324
                Tst = reaction_dict["Tst"]
                logKr = phase["logKr"]["values"][1]
                dG = R * Tst * log(10) * logKr
                reaction_dict["drsm_gibbs_energy"] = Dict{String, Any}(
                    "values" => [dG],
                    "units" => ["J/mol"]
                )

                reaction_dict["drsm_heat_capacity_p"] = Dict{String, Any}(
                    "values" => [""],
                    "units" => ["J/(mol*K)"]
                )

                reaction_dict["drsm_enthalpy"] = Dict{String, Any}(
                    "values" => [""],
                    "units" => ["J/mol"]
                )

                reaction_dict["drsm_entropy"] = Dict{String, Any}(
                    "values" => [""],
                    "units" => ["J/(mol*K)"]
                )

                # The Vm field in Phreeqc .dat is the molar volume of the solid phase in cm3/mol
                # reaction_dict["drsm_volume"] = Dict{String, Any}(
                #     "values" => [phase["drsm_volume"]],
                #     "units" => ["cm3/mol"]
                # )
                reaction_dict["drsm_volume"] = Dict{String, Any}(
                    "values" => [""],
                    "units" => ["J/bar"]
                )

                reaction_dict["datasources"] = ["Cemdata18"]

                push!(new_reactions_list, reaction_dict)
                println("New reaction added: $name")
            else
                @warn "Phase $name is missing required fields and will be skipped."
            end
        else
            println("Phase $name already exists in JSON, skipping.")
        end
    end

    append!(json_data["reactions"], new_reactions_list)
    return json_data
end

"""
    write_reaction(f, reaction)

Write a reaction Dict as JSON to a file (used for custom JSON output).

# Arguments
- `f`: An open IO stream to write to.
- `reaction`: The reaction dictionary to write.

# Returns
- Nothing. Writes directly to the file.
"""
function write_reaction(f, reaction)
    # Write the reaction as a JSON object to the file stream `f`
    write(f, "    {\n")
    write(f, "      \"symbol\": \"$(reaction["symbol"])\",\n")
    write(f, "      \"equation\": \"$(reaction["equation"])\",\n")
    if haskey(reaction, "comment")
        write(f, "      \"comment\": \"$(reaction["comment"])\",\n")
    end
    write(f, "      \"reactants\": [\n")
    for (j, reactant) in enumerate(reaction["reactants"])
        write(f, "        {\n")
        write(f, "          \"symbol\": \"$(reactant["symbol"])\",\n")
        write(f, "          \"coefficient\": $(reactant["coefficient"])\n")
        write(f, "        }")
        if j < length(reaction["reactants"])
            write(f, ",")
        end
        write(f, "\n")
    end
    write(f, "      ],\n")

    write(f, "      \"limitsTP\": {\n")
    write(f, "        \"range\": false,\n")
    write(f, "        \"lowerP\": 0.1,\n")
    write(f, "        \"lowerT\": 273.15,\n")
    write(f, "        \"upperP\": 1000000,\n")
    write(f, "        \"upperT\": 298.15\n")
    write(f, "      },\n")

    write(f, "      \"Tst\": 298.15,\n")
    write(f, "      \"Pst\": 100000,\n")

    write(f, "      \"TPMethods\": [\n")
    for (j, method) in enumerate(reaction["TPMethods"])
        write(f, "        {\n")
        if haskey(method, "method") && haskey(method["method"], "0")
            write(f, "          \"method\": {\n")
            write(f, "            \"0\": \"logk_fpt_function\"\n")
            write(f, "          },\n")
            write(f, "          \"limitsTP\": {\n")
            write(f, "            \"lowerP\": 0,\n")
            write(f, "            \"lowerT\": 273.15,\n")
            write(f, "            \"upperP\": 0,\n")
            write(f, "            \"upperT\": 273.15\n")
            write(f, "          },\n")
            write(f, "          \"logk_ft_coeffs\": {\n")
            write(f, "            \"values\": [\n")
            for (k, value) in enumerate(method["logk_ft_coeffs"]["values"])
                write(f, "              $(value)")
                if k < length(method["logk_ft_coeffs"]["values"])
                    write(f, ",")
                end
                write(f, "\n")
            end
            write(f, "            ]\n")
            write(f, "          }\n")
        elseif haskey(method, "method") && haskey(method["method"], "7")
            write(f, "          \"method\": {\n")
            write(f, "            \"7\": \"logk_3_term_extrap\"\n")
            write(f, "          }\n")
        elseif haskey(method, "method") && haskey(method["method"], "13")
            write(f, "          \"method\": {\n")
            write(f, "            \"13\": \"dr_volume_constant\"\n")
            write(f, "          }\n")
        end
        write(f, "        }")
        if j < length(reaction["TPMethods"])
            write(f, ",")
        end
        write(f, "\n")
    end
    write(f, "      ],\n")

    write(f, "      \"logKr\": {\n")
    write(f, "        \"values\": [$(reaction["logKr"]["values"][1])],\n")
    write(f, "        \"errors\": [2]\n")
    write(f, "      },\n")

    write(f, "      \"drsm_heat_capacity_p\": {\n")
    write(f, "        \"values\": [$(reaction["drsm_heat_capacity_p"]["values"][1])],\n")
    write(f, "        \"units\": [\"$(reaction["drsm_heat_capacity_p"]["values"][1])\"]\n")
    write(f, "      },\n")

    write(f, "      \"drsm_gibbs_energy\": {\n")
    write(f, "        \"values\": [$(reaction["drsm_gibbs_energy"]["values"][1])],\n")
    write(f, "        \"units\": [\"$(reaction["drsm_gibbs_energy"]["units"][1])\"]\n")
    write(f, "      },\n")

    write(f, "      \"drsm_enthalpy\": {\n")
    write(f, "        \"values\": [$(reaction["drsm_enthalpy"]["values"][1])],\n")
    write(f, "        \"units\": [\"$(reaction["drsm_enthalpy"]["units"][1])\"]\n")
    write(f, "      },\n")

    write(f, "      \"drsm_entropy\": {\n")
    write(f, "        \"values\": [$(reaction["drsm_entropy"]["values"][1])],\n")
    write(f, "        \"units\": [\"$(reaction["drsm_entropy"]["units"][1])\"]\n")
    write(f, "      },\n")

    write(f, "      \"drsm_volume\": {\n")
    write(f, "        \"values\": [$(reaction["drsm_volume"]["values"][1])],\n")
    write(f, "        \"units\": [\"$(reaction["drsm_volume"]["units"][1])\"]\n")
    write(f, "      },\n")

    write(f, "      \"datasources\": [\"Cemdata18\"]\n")

    write(f, "    }")
end

"""
    merge_json(json_path, dat_path, output_path)

Merge a Cemdata .dat file into a ThermoFun JSON, preserving field order.

# Arguments
- `json_path`: Path to the original JSON file.
- `dat_path`: Path to the Cemdata .dat file.
- `output_path`: Path to write the merged JSON.

# Returns
- Nothing. Writes the merged JSON to `output_path`.
"""
function merge_json(json_path, dat_path, output_path)
    # Read the initial JSON file
    initial_content = read(json_path, String)

    # Parse the initial JSON file to get field order
    json_data = JSON.parsefile(json_path)

    # Preserve the initial structure
    aqueous_species = get_aqueous_species(json_data)
    dat_content = read(dat_path, String)
    new_reactions = parse_phases(dat_content, aqueous_species)

    # Add new reactions
    merged_data = merge_reactions(json_data, new_reactions)

    # Write the output JSON file, preserving the initial order
    open(output_path, "w") do f
        # Find the start and end indices of the "reactions" section
        lines = split(initial_content, '\n')
        reactions_start = 0
        reactions_end = 0
        for (i, line) in enumerate(lines)
            if occursin("\"reactions\": [", line)
                reactions_start = i
            elseif reactions_start != 0 && occursin("\"elements\": [", line)
                reactions_end = i-1
                break
            end
        end

        # Write the initial content up to the start of reactions
        for i in 1:reactions_start-1
            write(f, lines[i] * "\n")
        end

        # Write the start line of reactions
        write(f, lines[reactions_start] * "\n")

        # Write all reactions (existing and new)
        for (i, reaction) in enumerate(merged_data["reactions"])
            write_reaction(f, reaction)
            if i < length(merged_data["reactions"])
                write(f, ",\n")
            else
                write(f, "\n")
            end
        end

        # Write the end of reactions and the rest of the file
        for i in reactions_end:length(lines)
            write(f, lines[i] * "\n")
        end
    end
end

"""
    read_thermofun(filename)

Read a ThermoFun JSON file and return DataFrames for elements, substances, and reactions.

# Arguments
- `filename`: Path to the ThermoFun JSON file.

# Returns
- `df_elements`: DataFrame of elements.
- `df_substances`: DataFrame of substances.
- `df_reactions`: DataFrame of reactions.
"""
function read_thermofun(filename; with_units=true, debug=false, all_properties=false)
    # 1. Read the JSON file
    data = open(JSON3.read, filename)

    function extract_scal_or_vect(x)
        if x === missing
            return x
        end
        if length(x) == 1
            return x[1]
        else
            return x
        end
    end

    # 2. Function to extract values, units, and errors from a field
    function extract_field_with_units(fields, field_name)
        if haskey(fields, field_name)
            values = get(fields[field_name], "values", missing)
            units = get(fields[field_name], "units", missing)
            errors = get(fields[field_name], "errors", missing)
            if !isa(values, Missing) && length(values) > 0
                return (values=extract_scal_or_vect(values), units=extract_scal_or_vect(units), errors=extract_scal_or_vect(errors))
            else
                return (values=missing, units=missing, errors=missing)
            end
        else
            return (values=missing, units=missing, errors=missing)
        end
    end

    # Function to extract coefficients and their metadata
    function extract_coeffs_with_units(method, coeff_name)
        if haskey(method, coeff_name)
            values = get(method[coeff_name], "values", missing)
            units = get(method[coeff_name], "units", missing)
            names = get(method[coeff_name], "names", missing)
            if !isa(values, Missing) && length(values) > 0
                return (values=extract_scal_or_vect(values), units=extract_scal_or_vect(units), names=extract_scal_or_vect(names))
            end
        end
        return (values=missing, units=missing, names=missing)
    end

    # Function to extract limitsTP as a NamedTuple
    function extract_limitsTP(limitsTP)
        if limitsTP === missing
            return (lowerT=missing, upperT=missing, lowerP=missing, upperP=missing, range=missing)
        else
            return (
                lowerT=get(limitsTP, "lowerT", missing),
                upperT=get(limitsTP, "upperT", missing),
                lowerP=get(limitsTP, "lowerP", missing),
                upperP=get(limitsTP, "upperP", missing),
                range=get(limitsTP, "range", missing)
            )
        end
    end

    # Function to parse a TPMethod object into a Dict
    function parse_TPMethod(method)
        method_dict = Dict()
        # Extract the method type
        method_dict["method"] = only(values(get(method, "method", Dict("" => missing))))
        # Extract limitsTP if present
        if haskey(method, "limitsTP")
            method_dict["limitsTP"] = extract_limitsTP(method["limitsTP"])
        end
        # Extract coefficients (e.g., logk_ft_coeffs) if present
        if haskey(method, "logk_ft_coeffs")
            method_dict["logk_ft_coeffs"] = extract_coeffs_with_units(method, "logk_ft_coeffs")
        end
        if haskey(method, "m_heat_capacity_ft_coeffs")
            method_dict["m_heat_capacity_ft_coeffs"] = extract_coeffs_with_units(method, "m_heat_capacity_ft_coeffs")
        end
        # Add other fields as needed
        return method_dict
    end

    # 3. Create the DataFrame for substances
    substances = data.substances
    df_substances = DataFrame(
        name = [get(s, "name", missing) for s in substances],
        symbol = [get(s, "symbol", missing) for s in substances],
        formula = [get(s, "formula", missing) for s in substances],
        # formula = [replace(get(s, "formula", missing), r"\|\-?\d+\|" => "") for s in substances],
        charge = [get(s, "formula_charge", missing) for s in substances],
        aggregate_state = [only(values(get(s, "aggregate_state", Dict("" => missing)))) for s in substances],
        class = [only(values(get(s, "class_", Dict("" => missing)))) for s in substances],
        Tst = [get(s, "Tst", missing) for s in substances],
        Pst = [get(s, "Pst", missing) for s in substances],
        TPMethods = [
            [parse_TPMethod(m) for m in get(s, "TPMethods", [])]
            for s in substances
        ],
        Cp = [extract_field_with_units(s, "sm_heat_capacity_p") for s in substances],
        ΔfG = [extract_field_with_units(s, "sm_gibbs_energy") for s in substances],
        ΔfH = [extract_field_with_units(s, "sm_enthalpy") for s in substances],
        S = [extract_field_with_units(s, "sm_entropy_abs") for s in substances],
        Vm = [extract_field_with_units(s, "sm_volume") for s in substances],
        datasources = [get(s, "datasources", missing) for s in substances],
    )
    complete_species_database!(df_substances; with_units=with_units, debug=debug, all_properties=all_properties)

    # 4. Create the DataFrame for reactions
    reactions = data.reactions
    df_reactions = DataFrame(
        symbol = [get(r, "symbol", missing) for r in reactions],
        equation = [get(r, "equation", missing) for r in reactions],
        reactants = [
            Dict(r["symbol"] => r["coefficient"] for r in get(r, "reactants", []))
            for r in reactions
        ],
        limitsTP = [extract_limitsTP(get(r, "limitsTP", missing)) for r in reactions],
        Tst = [get(r, "Tst", missing) for r in reactions],
        Pst = [get(r, "Pst", missing) for r in reactions],
        TPMethods = [
            [parse_TPMethod(m) for m in get(r, "TPMethods", [])]
            for r in reactions
        ],
        logKr = [extract_field_with_units(r, "logKr") for r in reactions],
        ΔrCp = [extract_field_with_units(r, "drsm_heat_capacity_p") for r in reactions],
        ΔrG = [extract_field_with_units(r, "drsm_gibbs_energy") for r in reactions],
        ΔrH = [extract_field_with_units(r, "drsm_enthalpy") for r in reactions],
        ΔrS = [extract_field_with_units(r, "drsm_entropy") for r in reactions],
        ΔrV = [extract_field_with_units(r, "drsm_volume") for r in reactions],
        datasources = [get(r, "datasources", missing) for r in reactions]
    )

    # 5. Create the DataFrame for elements
    elements = data.elements
    df_elements = DataFrame(
        symbol = [get(e, "symbol", missing) for e in elements],
        class = [only(values(get(e, "class_", Dict("" => missing)))) for e in elements],
        S = [extract_field_with_units(e, "entropy") for e in elements],
        atomic_mass = [extract_field_with_units(e, "atomic_mass") for e in elements],
        datasources = [get(e, "datasources", missing) for e in elements]
    )

    return df_elements, df_substances, df_reactions
end

function complete_species_database!(df_substances; with_units=true, debug=false, all_properties=false)
    colspecies = [try Species(s.formula;
                       name = s.name,
                       symbol = s.symbol,
                       aggregate_state = try eval(Meta.parse(s.aggregate_state)) catch; AS_UNDEF end,
                       class = try eval(Meta.parse(s.class)) catch; SC_UNDEF end,
                       )
                    catch; missing end for s in eachrow(df_substances)]

    insertcols!(df_substances, 1, :species => colspecies)

    if debug print_title("Property completion"; crayon=Crayon(foreground=:blue), style=:box, indent="") end

    function get_value(row, field::Symbol; debug=debug, crayon=Crayon(), with_units=with_units, default_unit=unit(1))            
        val = row[field].values
        vunit = row[field].units
        if debug>1 && (iszero(val) || (with_units && ismissing(vunit)))
            println(crayon("$(row.symbol) => $field=$val $vunit"))
            print(crayon"reset")
        end
        if with_units
            val = val*(try uparse(vunit) catch; default_unit end)
        end
        return val
    end

    function get_Cp_coef(row; debug=debug, crayon=Crayon(), with_units=with_units)            
        TPMethods = row.TPMethods
        idx = findfirst(d -> haskey(d, "method") && d["method"] == "cp_ft_equation", TPMethods)
        if !isnothing(idx)
            d = TPMethods[idx]
            tuple_coefs = d["m_heat_capacity_ft_coeffs"]
            values = tuple_coefs.values
            if debug>1 && !iszero(max(abs.(values[5:end])...)) println(crayon("$(row.symbol) => Cp=$values")) end
            if with_units
                units = tuple_coefs.units
                values = [values[i]*uparse(units[i]) for i=1:min(length(values), length(units))]
            end
            return values
        else
            return [get_value(row, :Cp; debug=debug, crayon=crayon, with_units=with_units, default_unit=J/(mol*K))]
        end
    end

    iters = debug ? eachrow(df_substances) : ProgressBar(eachrow(df_substances))

    for row in iters
        s = row.species
        if debug println(s) end
        Tref = with_units ? row.Tst*K : row.Tst
        s.Tref = Tref
        if !with_units
            s.molar_mass = ustrip(s.molar_mass)
        end

        if all_properties
            coeffa = float.(get_Cp_coef(row; debug=debug, crayon=crayon"green"))
            s.Cp = ThermoFunction(:Cp, coeffa; Tref=Tref)

            ΔfH0 = get_value(row, :ΔfH; debug=debug, crayon=crayon"red", with_units=with_units, default_unit=J/mol)
            coeffaH = [float(zero(ΔfH0)); coeffa]
            fH = ThermoFunction(:H, coeffaH; Tref=Tref)
            coeffaH[1] = ΔfH0-Base.invokelatest(fH, Tref)
            s.ΔfH = ThermoFunction(:H, coeffaH; Tref=Tref)

            S0 = get_value(row, :S; debug=debug, crayon=crayon"red", with_units=with_units, default_unit=J/(mol*K))
            coeffaS = [float(zero(S0)); coeffa]
            fS = ThermoFunction(:S, coeffaS; Tref=Tref)
            coeffaS[1] = S0-Base.invokelatest(fS, Tref)
            s.S = ThermoFunction(:S, coeffaS; Tref=Tref)

            ΔfG0 = get_value(row, :ΔfG; debug=debug, crayon=crayon"blue", with_units=with_units, default_unit=J/mol)
            coeffaG = [float(zero(ΔfG0)); float(zero(S0)); -coeffa]
            fG = ThermoFunction(:G, coeffaG; Tref=Tref)
            coeffaG[1] = ΔfG0+coeffaS[1]*Tref-Base.invokelatest(fG, Tref)
            coeffaG[2] = -coeffaS[1]
            s.ΔfG = ThermoFunction(:G, coeffaG; Tref=Tref)

            Vm = uconvert(cm^3 , get_value(row, :Vm; crayon=crayon"blue", with_units=true, default_unit=J/bar))
            s.Vm = with_units ? Vm/mol : ustrip(Vm)
        end
    end

    colcemspecies = [try CemSpecies(s) catch; missing end for s in colspecies]                
    insertcols!(df_substances, 2, :cemspecies => colcemspecies)

    return df_substances

end

function build_species_database(df_substances; kwargs...)
    return complete_species_database!(copy(df_substances); kwargs...)
end


"""
    extract_primary_species(file_path)

Extract primary species from a Cemdata .dat file, returning a DataFrame.

# Arguments
- `file_path`: Path to the Cemdata .dat file.

# Returns
- DataFrame of primary species.
"""
function extract_primary_species(file_path)
    # Lire le fichier
    lines = readlines(file_path)

    # Trouver les indices des sections
    start_idx = 0
    end_idx = 0
    in_primary_section = false

    for (i, line) in enumerate(lines)
        stripped = strip(line)

        # Détecter le début de la section
        if startswith(stripped, "SOLUTION_SPECIES")
            start_idx = i + 1  # Commencer après SOLUTION_SPECIES
            in_primary_section = true

            # Trouver la première espèce (ligne avec "=")
            while start_idx <= length(lines)
                next_line = strip(lines[start_idx])
                if startswith(next_line, "# PMATCH MASTER SPECIES") ||
                   occursin("=", next_line)
                    break
                end
                start_idx += 1
            end

            # Si on a trouvé "# PMATCH MASTER SPECIES", commencer après
            if startswith(strip(lines[start_idx]), "# PMATCH MASTER SPECIES")
                start_idx += 1
            end
        end

        # Détecter la fin de la section
        if in_primary_section && startswith(stripped, "# PMATCH SECONDARY MASTER SPECIES")
            end_idx = i - 1
            break
        end
    end

    # Si on n'a pas trouvé la fin, prendre jusqu'à la fin du fichier
    if in_primary_section && end_idx == 0
        end_idx = length(lines)
    end

    # Extraire les données
    species_data = []

    for i in start_idx:end_idx
        line = strip(lines[i])

        # Traiter les espèces (lignes avec "=")
        if occursin("=", line)
            parts = split(line, "=")
            current_species = strip(parts[1])

            # Conversion au format ThermoFun
            if current_species == "e-"
                symbol = "Zz"
            elseif occursin(r"[\+\-]\d*$", current_species)
                symbol = current_species
            else
                symbol = current_species * "@"
            end

            push!(species_data, (species=current_species, symbol=symbol, gamma=Float64[]))
        end

        # Traiter les paramètres gamma (lignes avec "-gamma")
        if startswith(line, "-gamma") && !isempty(species_data)
            parts = split(line)
            gamma_values = Float64[]
            for val in parts[2:end]  # Commencer après "-gamma"
                num = tryparse(Float64, val)
                if num !== nothing
                    push!(gamma_values, num)
                end
            end

            if !isempty(gamma_values)
                # Mettre à jour la dernière espèce
                last_entry = species_data[end]
                species_data[end] = (species=last_entry.species,
                                    symbol=last_entry.symbol,
                                    gamma=gamma_values)
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

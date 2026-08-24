# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using JSON
using OrderedCollections

"""
    merge_reactions(json_data::Dict, new_reactions::Dict) -> Dict

Merge new reactions from PHREEQC .dat into existing ThermoFun JSON data.

# Arguments

  - `json_data`: existing ThermoFun database as a dictionary.
  - `new_reactions`: dictionary of new phase reactions to add.

# Returns

  - Updated `json_data` with new reactions appended.

Only adds reactions that don't already exist (by symbol) and that have complete
required fields (logKr, analytical_expression, equation).
"""
function merge_reactions(json_data, new_reactions)
    existing_symbols = Set{String}()
    for reaction in json_data["reactions"]
        push!(existing_symbols, reaction["symbol"])
    end

    new_reactions_list = []
    for (name, phase) in new_reactions
        if !(name in existing_symbols)
            if haskey(phase, "logKr") &&
                    haskey(phase, "analytical_expression") &&
                    haskey(phase, "equation")
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
                    "upperT" => 298.15,
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
                            "upperT" => 273.15,
                        ),
                        "logk_ft_coeffs" => Dict{String, Any}(
                            "values" => vcat(
                                phase["analytical_expression"],
                                zeros(max(0, 12 - length(phase["analytical_expression"]))),
                            ),
                        ),
                    ),
                    Dict{String, Any}("method" => Dict("7" => "logk_3_term_extrap")),
                    Dict{String, Any}("method" => Dict("13" => "dr_volume_constant")),
                ]

                reaction_dict["logKr"] = phase["logKr"]

                R = 8.31446261815324
                Tst = reaction_dict["Tst"]
                logKr = phase["logKr"]["values"][1]
                dG = R * Tst * log(10) * logKr
                reaction_dict["drsm_gibbs_energy"] = Dict{String, Any}(
                    "values" => [dG], "units" => ["J/mol"]
                )

                reaction_dict["drsm_heat_capacity_p"] = Dict{String, Any}(
                    "values" => [""], "units" => ["J/(mol*K)"]
                )

                reaction_dict["drsm_enthalpy"] = Dict{String, Any}(
                    "values" => [""], "units" => ["J/mol"]
                )

                reaction_dict["drsm_entropy"] = Dict{String, Any}(
                    "values" => [""], "units" => ["J/(mol*K)"]
                )

                reaction_dict["drsm_volume"] = Dict{String, Any}(
                    "values" => [""], "units" => ["J/bar"]
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
    write_reaction(f::IO, reaction::Dict)

Write a single reaction dictionary to an IO stream in JSON format.

# Arguments

  - `f`: output IO stream.
  - `reaction`: reaction dictionary with all required ThermoFun fields.

Helper function used by `merge_json` to write formatted JSON output.
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
    write(f, "        \"units\": [\"$(reaction["drsm_heat_capacity_p"]["units"][1])\"]\n")
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

    return write(f, "    }")
end

"""
    merge_json(json_path::AbstractString, dat_path::AbstractString, output_path::AbstractString)

Merge PHREEQC .dat phase data into a ThermoFun JSON database file.

# Arguments

  - `json_path`: path to input ThermoFun JSON file.
  - `dat_path`: path to PHREEQC .dat file containing phase definitions.
  - `output_path`: path for output merged JSON file.

Reads both files, extracts phases from the .dat file, merges them into the JSON
database structure, and writes the result preserving the original JSON formatting.
"""
function merge_json(json_path, dat_path, output_path)
    # The two inputs are resolved against the bundled `data/` directory when the
    # working directory does not hold them; `output_path` never is, because it
    # names a file to create.
    json_path = resolve_data_path(json_path)
    dat_path = resolve_data_path(dat_path)

    # Read the initial JSON file
    initial_content = read(json_path, String)

    # Parse the initial JSON file to get field order
    json_data = JSON.parsefile(json_path)

    # Preserve the initial structure
    dat_content = read(dat_path, String)
    new_reactions = parse_phases(dat_content)

    # Add new reactions
    merged_data = merge_reactions(json_data, new_reactions)

    # Write the output JSON file, preserving the initial order
    return open(output_path, "w") do f
        # Find the start and end indices of the "reactions" section
        lines = split(initial_content, '\n')
        reactions_start = 0
        reactions_end = 0
        for (i, line) in enumerate(lines)
            if occursin("\"reactions\": [", line)
                reactions_start = i
            elseif reactions_start != 0 && occursin("\"elements\": [", line)
                reactions_end = i - 1
                break
            end
        end

        # Write the initial content up to the start of reactions
        for i in 1:(reactions_start - 1)
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

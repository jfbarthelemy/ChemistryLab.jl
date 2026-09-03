# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using OrderedCollections
using Unicode

"""
    stoich_coef_round(x::T; tol=1e-4) where {T<:Real} -> Union{Int, Rational, Float64}
    stoich_coef_round(x) -> Any

Round stoichiometric coefficients to integer, rational, or float representation.

# Arguments

  - `x`: numeric value to round.
  - `tol`: tolerance for rounding decisions (default 1e-4).

# Returns

  - Integer if close to a whole number, Rational if a simple fraction (denominator < 10),
    or Float64 rounded to 5 digits otherwise. Non-numeric inputs are returned unchanged.

# Examples

```jldoctest
julia> stoich_coef_round(2.00001)
2

julia> stoich_coef_round(0.3333)
1//3

julia> stoich_coef_round(3.14159)
3.14159
```
"""
function stoich_coef_round(x::T; tol = 1.0e-3) where {T <: Real}
    try
        if isapprox(x, round(x); atol = tol)
            return Int(round(x))
        end

        rat = rationalize(x; tol = tol)
        if isapprox(x, float(rat); atol = tol)
            if 1 < denominator(rat) < 10
                return rat
            end
        end

        return round(x; digits = 5)
    catch e
        return x
    end
end

stoich_coef_round(x) = x

"""
    phreeqc_to_unicode(s::AbstractString) -> String

Convert a PHREEQC formula string to Unicode representation with subscripts and superscripts.

# Arguments

  - `s`: PHREEQC-formatted chemical formula string.

# Returns

  - Unicode-formatted string with subscript coefficients and superscript charges.

# Examples

```jldoctest
julia> phreeqc_to_unicode("Ca+2")
"Ca²⁺"

julia> phreeqc_to_unicode("SO4-2")
"SO₄²⁻"
```
"""
function phreeqc_to_unicode(s::AbstractString)
    chars = collect(s)

    ind_sign = findall(
        i -> (chars[i] == '+' || chars[i] == '-') && (i > 1 && chars[i - 1] != ' '),
        1:length(chars),
    )

    for i in ind_sign
        sign = chars[i]
        j = i
        while j < length(chars) && (isnumeric(chars[j + 1]) || chars[j + 1] == '.')
            chars[j] = dict_normal_to_super[chars[j + 1]]
            j += 1
        end
        chars[j] = dict_normal_to_super[sign]
    end

    s = join(chars)

    s = replace(s, r"-?\d+\.?\d*" => x -> string(stoich_coef_round(parse(Float64, x))))

    matches = collect(eachmatch(r"(\d+)\/\/(\d+)", s))
    for m in reverse(matches)
        num = tryparse(Int, m.captures[1])
        den = tryparse(Int, m.captures[2])
        rat = stoich_coef_round(num // den)
        replacement = get(dict_frac_unicode, rat, string(rat))
        start_idx = m.offset
        end_idx = start_idx + sizeof(m.match) - 1
        prefix = if start_idx > first(eachindex(s))
            s[first(eachindex(s)):prevind(s, start_idx)]
        else
            ""
        end
        suffix = end_idx < sizeof(s) ? s[(end_idx + 1):end] : ""
        s = prefix * replacement * suffix
    end

    chars = collect(s)

    ind_sign = findall(
        i ->
        chars[i] in keys(dict_normal_to_sub) &&
            i > 1 &&
            chars[i - 1] != ' ' &&
            !(chars[i - 1] in keys(dict_normal_to_sub)),
        1:length(chars),
    )

    for i in ind_sign
        j = i
        while j <= length(chars) && chars[j] in keys(dict_normal_to_sub)
            chars[j] = dict_normal_to_sub[chars[j]]
            j += 1
        end
    end

    return join(chars)
end

"""
    unicode_to_phreeqc(s::AbstractString) -> String

Convert a Unicode formula string back to PHREEQC format.

# Arguments

  - `s`: Unicode-formatted chemical formula.

# Returns

  - PHREEQC-formatted string with plain text charges and coefficients.

# Examples

```jldoctest
julia> unicode_to_phreeqc("Ca²⁺")
"Ca+2"

julia> unicode_to_phreeqc("SO₄²⁻")
"SO4-2"
```
"""
function unicode_to_phreeqc(s::AbstractString)
    chars = collect(s)

    ind_sign = findall(k -> k == '⁺' || k == '⁻', chars)
    for i in ind_sign
        sign = chars[i]
        j = i
        while j > 1 && chars[j - 1] in keys(dict_super_to_normal)
            chars[j] = chars[j - 1]
            j -= 1
        end
        chars[j] = sign
    end
    s = super_to_normal(join(chars))

    s = replace(s, r"[₀₁₂₃₄₅₆₇₈₉]" => c -> dict_sub_to_normal[c[1]])

    pattern = r"(\d+(\.\d+)?|)([¼½¾⅐⅑⅒⅓⅔⅕⅖⅗⅘⅙⅚⅛⅜⅝⅞])"
    matches = collect(eachmatch(pattern, s))
    for m in reverse(matches)
        float_part = isempty(m.captures[1]) ? 0.0 : parse(Float64, m.captures[1])
        sum_value = float_part + dict_unicode_frac[m.captures[3][1]]
        replacement = string(stoich_coef_round(sum_value))
        start_idx = m.offset
        end_idx = start_idx + sizeof(m.match) - 1
        prefix = if start_idx > first(eachindex(s))
            s[first(eachindex(s)):prevind(s, start_idx)]
        else
            ""
        end
        suffix = end_idx < sizeof(s) ? s[(end_idx + 1):end] : ""
        s = prefix * replacement * suffix
    end

    return s
end

"""
    merge_upper_lower(graphemes::Vector{<:AbstractString}) -> Vector{String}

Merge consecutive graphemes where an uppercase letter is followed by a lowercase letter.

# Arguments

  - `graphemes`: vector of grapheme strings.

# Returns

  - Vector with merged element symbols (e.g., ["C", "a"] becomes ["Ca"]).
"""
function merge_upper_lower(graphemes::Vector{<:AbstractString})
    result = String[]
    i = 1
    while i <= length(graphemes)
        current = graphemes[i]
        if i < length(graphemes)
            last_char = current[end]
            next_first_char = graphemes[i + 1][1]
            if isuppercase(last_char) && islowercase(next_first_char)
                current *= graphemes[i + 1]
                i += 1
            end
        end
        push!(result, current)
        i += 1
    end
    return result
end

"""
    colored_formula(s::AbstractString; colorcharge=true) -> String

Generate a terminal-colored representation of a chemical formula.

# Arguments

  - `s`: formula string.
  - `colorcharge`: if true, color the charge portion (default true).

# Returns

  - String with ANSI color codes for terminal display.
"""
function colored_formula(s::AbstractString; colorcharge = true)
    superscript_digits = ["⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹", "⁺", "⁻"]

    colored_graph = merge_upper_lower(collect(graphemes(s)))
    ind_sign = findlast(k -> k == "⁺" || k == "⁻" || k == "+" || k == "-", colored_graph)
    if isnothing(ind_sign)
        ind_sign = length(colored_graph) + 1
    end

    idx_sign = Int[]
    idx_atoms = Int[]
    idx_par = Int[]
    for (i, c) in enumerate(colored_graph)
        if colorcharge && (i >= ind_sign || c in superscript_digits)
            push!(idx_sign, i)
        end
        if Symbol(c) in ATOMIC_ORDER || Symbol(c) in OXIDE_ORDER
            push!(idx_atoms, i)
        end
        if c in ["(", ")", "[", "]", "{", "}", "@", "|"]
            push!(idx_par, i)
        end
    end
    idx_stoich = setdiff(1:length(colored_graph), union(idx_sign, idx_atoms, idx_par))
    colored_graph[idx_sign] .= string.(COL_CHARGE.(colored_graph[idx_sign]))
    colored_graph[idx_par] .= string.(COL_PAR.(colored_graph[idx_par]))
    colored_graph[idx_stoich] .= string.(COL_STOICH_INT.(colored_graph[idx_stoich]))

    return join(colored_graph)
end

"""
    parse_formula(formula::AbstractString) -> OrderedDict{Symbol,Number}

Parse a chemical formula string into an atomic composition dictionary.

# Arguments

  - `formula`: formula string (supports parentheses, brackets, rational/decimal coefficients).

# Returns

  - OrderedDict mapping element symbols to stoichiometric coefficients.

# Examples

```jldoctest
julia> parse_formula("H2O")
OrderedDict{Symbol, Int64} with 2 entries:
  :H => 2
  :O => 1

julia> parse_formula("Ca(OH)2")
OrderedDict{Symbol, Int64} with 3 entries:
  :Ca => 1
  :O  => 2
  :H  => 2
```
"""
function parse_formula(formula::AbstractString)
    function safe_nextind(s::AbstractString, i::Integer, n::Integer = 1)
        last_i = lastindex(s)
        idx = i
        for _ in 1:n
            if idx > last_i
                return last_i + 1
            end
            idx = nextind(s, idx)
        end
        return idx
    end

    formula = replace(formula, ":" => "", "{" => "(", "}" => ")", "[" => "(", "]" => ")")
    formula = replace(formula, r"\|\-?\d+\|" => "")
    formula = replace(formula, r"\|" => "")
    formula = unicode_to_phreeqc(String(formula))

    counts = OrderedDict{Symbol, Number}()

    i = firstindex(formula)
    while i <= lastindex(formula)
        c = formula[i]
        if c == '('
            depth = 1
            j = safe_nextind(formula, i)
            while j <= lastindex(formula) && depth > 0
                if formula[j] == '('
                    depth += 1
                elseif formula[j] == ')'
                    depth -= 1
                end
                j = safe_nextind(formula, j)
            end
            inner = formula[safe_nextind(formula, i):prevind(formula, j)]

            rest = j <= lastindex(formula) ? formula[j:end] : ""

            m = match(r"^([0-9]+//[0-9]+|[0-9]+(?:\.[0-9]+)?)", rest)
            factor = if (m === nothing)
                1
            else
                begin
                    s = m.match
                    factor =
                        occursin("//", s) ? parse(Rational{Int}, s) : parse(Float64, s)
                end
            end
            offset = (m === nothing) ? 0 : length(m.match)

            for (el, n) in parse_formula(inner)
                counts[el] = get(counts, el, 0) + n * factor
            end

            i = safe_nextind(formula, j, offset)
        else
            m = match(
                r"^(\p{Lu}[\p{Ll}\u0300-\u036F]?)(([0-9]+//[0-9]+)|([0-9]+(?:\.[0-9]+)?))?",
                formula[i:end],
            )
            if m !== nothing
                el, countstr = m.captures
                el = Symbol(el)
                cnt = if countstr === nothing || isempty(countstr)
                    1
                elseif occursin("//", countstr)
                    parse(Rational{Int}, countstr)
                else
                    stoich_coef_round(parse(Float64, countstr))
                end

                if cnt isa Rational && denominator(cnt) == 1
                    cnt = Int(numerator(cnt))
                elseif cnt isa Float64 && isinteger(cnt)
                    cnt = Int(cnt)
                end

                counts[el] = get(counts, el, 0) + cnt

                i = safe_nextind(formula, i, length(m.match))
            else
                i = safe_nextind(formula, i)
            end
        end
    end

    return OrderedDict(k => stoich_coef_round(v) for (k, v) in counts)
end

"""
    extract_charge(formula::AbstractString) -> Int

Extract the formal charge from a chemical formula string.

# Arguments

  - `formula`: formula string with optional charge notation (e.g., "+2", "-", "3+").

# Returns

  - Integer charge value (0 if no charge present).

# Examples

```jldoctest
julia> extract_charge("Ca+2")
2

julia> extract_charge("SO4-2")
-2

julia> extract_charge("H2O")
0
```
"""
function extract_charge(formula::AbstractString)
    m = match(r"([+-])([0-9]*)$", unicode_to_phreeqc(formula))
    if m === nothing
        return 0
    else
        sign = m.captures[1] == "+" ? 1 : -1
        val = m.captures[2] == "" ? 1 : parse(Int, m.captures[2])
        return sign * val
    end
end

"""
    to_mendeleev(oxides::AbstractDict{Symbol,T}) where {T<:Number} -> OrderedDict{Symbol,Number}

Convert cement oxide notation to Mendeleev element composition.

# Arguments

  - `oxides`: dictionary mapping oxide symbols (C, S, A, etc.) to coefficients.

# Returns

  - OrderedDict mapping element symbols to stoichiometric coefficients.

# Examples

```jldoctest
julia> to_mendeleev(OrderedDict(:C => 1, :S => 2))
OrderedDict{Symbol, Int64} with 3 entries:
  :Ca => 1
  :O  => 5
  :Si => 2
```

# See also

  - [`CEMENT_TO_MENDELEEV`](@ref)
"""
function to_mendeleev(oxides::AbstractDict{Symbol, T}) where {T <: Number}
    result = OrderedDict{Symbol, Number}()
    for (ox, coef) in oxides
        if ox ∉ [:Zz, :Zz⁺, :e, :e⁻]
            idx = findfirst(p -> p.first == ox, CEMENT_TO_MENDELEEV)
            idx !== nothing || error("$(ox) is not a valid oxide identifier")
            mend = CEMENT_TO_MENDELEEV[idx].second
            for (k, v) in mend
                result[k] = get(result, k, 0) + v * coef
            end
        end
    end
    return if length(result) > 0
        OrderedDict(k => stoich_coef_round(v) for (k, v) in result)
    else
        result
    end
end

"""
    parse_equation(equation::AbstractString) -> Tuple{OrderedDict{String,Real}, OrderedDict{String,Real}, Char}

Parse a chemical equation string into reactants, products, and equality sign.

# Arguments

  - `equation`: equation string with reactants and products separated by an equality operator.

# Returns

  - Tuple of (reactants dict, products dict, equal sign char).

# Examples

```jldoctest
julia> reactants, products, sign = parse_equation("2H2 + O2 = 2H2O")
(OrderedDict("H2" => 2, "O2" => 1), OrderedDict("H2O" => 2), '=')

julia> reactants["H2"]
2

julia> products["H2O"]
2
```
"""
function parse_equation(equation::AbstractString)
    equal_sign = '='
    for c in equation
        if c in EQUAL_REACTION_SET
            equal_sign = c
            break
        end
    end

    sides = strip.(split(equation, EQUAL_REACTION))
    nsides = length(sides)
    left_side = nsides > 0 ? sides[1] : ""
    right_side = nsides > 1 ? sides[2] : ""

    function parse_side(side::AbstractString)
        terms = split(side, " +")
        result = OrderedDict{String, Real}()

        for term in terms
            t = strip(term)

            m = match(r"^(?<coeff>[-+]?\d+//\d+|[-+]?\d*\.?\d+)?\s*(?<formula>.+)$", t)

            if m !== nothing
                coeff_str = m[:coeff]
                formula = strip(m[:formula])

                coeff = if coeff_str === nothing || coeff_str == ""
                    1
                else
                    eval(Meta.parse(coeff_str))
                end

                if !(coeff isa Real)
                    error("Invalid coefficient: $coeff_str")
                end

                result[formula] = coeff
            else
                error("Unexpected term format: $term")
            end
        end

        return OrderedDict(k => stoich_coef_round(v) for (k, v) in result)
    end

    reactants = if left_side == "∅" || left_side == ""
        OrderedDict{String, Int}()
    else
        parse_side(left_side)
    end
    products = if right_side == "∅" || right_side == ""
        OrderedDict{String, Int}()
    else
        parse_side(right_side)
    end

    return reactants, products, equal_sign
end

"""
    colored_equation(equation::AbstractString) -> String

Generate a terminal-colored representation of a chemical equation.

# Arguments

  - `equation`: chemical equation string.

# Returns

  - String with ANSI color codes for terminal display.
"""
function colored_equation(equation::AbstractString)
    reactants, products, equal_sign = parse_equation(equation)
    left_side = if isempty(reactants)
        "∅"
    else
        join(
            [
                string(
                    COL_STOICH_EXT(
                        if isone(v)
                            ""
                        elseif v < 0
                            "($(v))"
                        else
                            string(v)
                        end,
                    )
                ) * colored_formula(k) for (k, v) in reactants
            ],
            " + ",
        )
    end
    right_side = if isempty(products)
        "∅"
    else
        join(
            [
                string(
                    COL_STOICH_EXT(
                        if isone(v)
                            ""
                        elseif v < 0
                            "($(v))"
                        else
                            string(v)
                        end,
                    )
                ) * colored_formula(k) for (k, v) in products
            ],
            " + ",
        )
    end
    return left_side * " " * string(COL_PAR(string(equal_sign))) * " " * right_side
end

"""
    format_equation(coeffs::AbstractDict; scaling=1, equal_sign='=') -> String

Format a stoichiometric coefficient dictionary into an equation string.

# Arguments

  - `coeffs`: dictionary mapping species to stoichiometric coefficients (negative = reactants, positive = products).
  - `scaling`: optional scaling factor for all coefficients (default 1).
  - `equal_sign`: equality operator to use (default '=').

# Returns

  - Formatted equation string with automatic electron balancing if needed.

# Examples

```jldoctest
julia> coeffs = OrderedDict("H2" => -2, "O2" => -1, "H2O" => 2);

julia> format_equation(coeffs)
"2H2 + O2 = 2H2O"
```
"""
function format_equation(coeffs::AbstractDict; scaling = 1, equal_sign = '=')
    # Separate reactants and products
    reactants = String[]
    products = String[]
    total_charge_left = 0
    total_charge_right = 0

    for (species, coeff) in coeffs
        if species !== "Zz"
            coeff = stoich_coef_round(coeff * scaling)

            # Format the coefficient
            abs_coeff = coeff < 0 ? -coeff : coeff
            coeff_str = if isapprox(abs_coeff, 1; atol = 1.0e-6)
                ""
            elseif isinteger(abs_coeff)
                string(Int(abs_coeff))
            else
                string(abs_coeff)
            end

            if coeff > 0
                push!(products, "$coeff_str$species")
                total_charge_right += coeff * extract_charge(species)
            elseif coeff < 0
                push!(reactants, "$coeff_str$species")
                total_charge_left += coeff * extract_charge(species)
            elseif coeff == 0
            else
                push!(products, "($coeff_str)$species")
                total_charge_right += coeff * extract_charge(species)
            end
        end
    end

    # Build the initial equation
    left_side = join(reactants, " + ")
    right_side = join(products, " + ")

    # Compute the charge difference (corrected)
    charge_diff = total_charge_right + total_charge_left

    # Balance charges if necessary
    if !isapprox(charge_diff, 0; atol = 1.0e-6)
        needed_e = stoich_coef_round(abs(charge_diff))
        e_term = needed_e == 1 ? "e⁻" : "$needed_e" * "e⁻"

        if charge_diff < 0
            # Add e- to the left (reactants)
            left_side = isempty(left_side) ? e_term : "$left_side + $e_term"
        else
            # Add e- to the right (products)
            right_side = isempty(right_side) ? e_term : "$right_side + $e_term"
        end
    end

    if length(left_side) == 0
        left_side = "∅"
    end
    if length(right_side) == 0
        right_side = "∅"
    end

    return "$left_side $(isnothing(equal_sign) ? '=' : equal_sign) $right_side"
end

"""
    add_parentheses_if_needed(s::String) -> String

Add parentheses to an expression string if it contains root-level +/- operators.

# Arguments

  - `s`: expression string.

# Returns

  - String wrapped in parentheses if needed, unchanged otherwise.

Parentheses are added only when the expression contains addition or subtraction
operators at the root level (outside any existing parentheses).
"""
function add_parentheses_if_needed(s::String)
    # Return early if s is already parenthesized fully and balanced
    if startswith(s, "(") && endswith(s, ")")
        # Check if outer parentheses fully enclose s
        count = 0
        for (i, c) in enumerate(s)
            if c == '('
                count += 1
            elseif c == ')'
                count -= 1
                if count == 0 && i != lastindex(s)
                    break
                end
            end
        end
        if count == 0
            return s  # already properly parenthesized
        end
    end

    # Helper to detect + or - at root level (outside any parentheses)
    function has_root_level_plusminus(str)
        lvl = 0
        for ch in str
            if ch == '('
                lvl += 1
            elseif ch == ')'
                lvl -= 1
            elseif (ch == '+' || ch == '-') && lvl == 0
                return true
            end
        end
        return false
    end

    # Add parentheses only if necessary
    if has_root_level_plusminus(s)
        return "(" * s * ")"
    else
        return s
    end
end

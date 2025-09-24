const ATOMIC_ORDER = [
    :Ca, :Na, :K, :Mg, :Sr, :Ba,  # Alkaline earth/alkali metals
    :Al, :Si, :Fe, :Ti, :Mn, :Cr,  # Transition metals
    :S, :P, :C, :B,               # Non-metals
    :H, :O,                       # H before O for H₂O, OH⁻, etc.
    :F, :Cl, :Br, :I,             # Halogens
    :N,                           # Nitrogen
    :Zz                            # Charge
]

# const cement_to_mendeleev = [
#     "C" => "CaO",  
#     "S" => "SiO2", 
#     "A" => "Al2O3",
#     "F" => "Fe2O3",
#     "K" => "K2O",  
#     "N" => "Na2O", 
#     "M" => "MgO",  
#     "P" => "P2O5", 
#     "T" => "TiO2", 
#     "C̄" => "CO2",  
#     "S̄" => "SO3",  
#     "N̄" => "NO3",  
#     "H" => "H2O",  
# ]

const cement_to_mendeleev = [
    :C  => Dict(:Ca => 1, :O => 1),
    :S  => Dict(:Si => 1, :O => 2),
    :A  => Dict(:Al => 2, :O => 3),
    :F  => Dict(:Fe => 2, :O => 3),
    :K  => Dict(:K  => 2, :O => 1),
    :N  => Dict(:Na => 2, :O => 1),
    :M  => Dict(:Mg => 1, :O => 1),
    :P  => Dict(:P  => 2, :O => 5),
    :T  => Dict(:Ti => 1, :O => 2),
    :C̄ => Dict(:C  => 1, :O => 2),
    :S̄ => Dict(:S  => 1, :O => 3),
    :N̄ => Dict(:N  => 1, :O => 3),
    :H  => Dict(:H  => 2, :O => 1),
]

const OXIDE_ORDER = [
    :C, :S, :A, :F, :K, :N, :M, :P, :T, :C̄, :S̄, :N̄, :H
]

const dict_super_to_normal = Dict{Char,Char}(
    '⁰' => '0',
    '¹' => '1',
    '²' => '2',
    '³' => '3',
    '⁴' => '4',
    '⁵' => '5',
    '⁶' => '6',
    '⁷' => '7',
    '⁸' => '8',
    '⁹' => '9',
    '⁺' => '+',
    '⁻' => '-',
    '.' => '.',
)

const dict_normal_to_super = Dict{Char,Char}(
    '0' => '⁰',
    '1' => '¹',
    '2' => '²',
    '3' => '³',
    '4' => '⁴',
    '5' => '⁵',
    '6' => '⁶',
    '7' => '⁷',
    '8' => '⁸',
    '9' => '⁹',
    '+' => '⁺',
    '-' => '⁻',
    '.' => '.',
)

const dict_sub_to_normal = Dict{Char,Char}(
    '₀' => '0',
    '₁' => '1',
    '₂' => '2',
    '₃' => '3',
    '₄' => '4',
    '₅' => '5',
    '₆' => '6',
    '₇' => '7',
    '₈' => '8',
    '₉' => '9',
    '.' => '.',
)

const dict_normal_to_sub = Dict{Char,Char}(
    '0' => '₀',
    '1' => '₁',
    '2' => '₂',
    '3' => '₃',
    '4' => '₄',
    '5' => '₅',
    '6' => '₆',
    '7' => '₇',
    '8' => '₈',
    '9' => '₉',
    '.' => '.',
)

const dict_frac_unicode = Dict(
        1//4 => "¼", 1//2 => "½", 3//4 => "¾",
        1//7 => "⅐", 1//9 => "⅑", 1//10 => "⅒",
        1//3 => "⅓", 2//3 => "⅔",
        1//5 => "⅕", 2//5 => "⅖", 3//5 => "⅗", 4//5 => "⅘",
        1//6 => "⅙", 5//6 => "⅚",
        1//8 => "⅛", 3//8 => "⅜", 5//8 => "⅝", 7//8 => "⅞"
    )

const dict_unicode_frac = Dict(
            '¼' => 1//4, '½' => 1//2, '¾' => 3//4,
            '⅐' => 1//7, '⅑' => 1//9, '⅒' => 1//10,
            '⅓' => 1//3, '⅔' => 2//3,
            '⅕' => 1//5, '⅖' => 2//5, '⅗' => 3//5, '⅘' => 4//5,
            '⅙' => 1//6, '⅚' => 5//6,
            '⅛' => 1//8, '⅜' => 3//8, '⅝' => 5//8, '⅞' => 7//8
        )

"Return whether `c` is a numeric superscript or ⁺/⁻."
issuperscript(c::Char) = c in keys(superscripts)

"Return whether `c` is a numeric subscript."
issubscript(c::Char) = c in keys(subscripts)

"Convert all numeric superscripts or ⁺/⁻ in `s` to normal line."
super_to_normal(s::AbstractString) = replace(s, dict_super_to_normal...)

"Convert all normal characters or +/- in `s` to numeric superscripts."
normal_to_super(s::AbstractString) = replace(s, dict_normal_to_super...)

"Convert all numeric subscripts in `s` to normal line."
sub_to_normal(s::AbstractString) = replace(s, dict_sub_to_normal...)

"Convert all normal characters in `s` to numeric subscripts ."
normal_to_sub(s::AbstractString) = replace(s, dict_normal_to_sub...)

function stoich_coef_round(x::T; tol=1e-4) where {T<:Real}
    try
        if isapprox(x, round(x); atol=tol)
            return Int(round(x))
        end

        rat = rationalize(x; tol=tol)
        if isapprox(x, float(rat); atol=tol)
            if 1<denominator(rat)<10
                return rat
            end
        end
        
        return round(x; digits=5)
    catch e
        return x
    end
end

stoich_coef_round(x) = x

function phreeqc_to_unicode(s::AbstractString)
    chars = collect(s)

    ind_sign = findall(i -> (chars[i] == '+' || chars[i] == '-') && (i > 1 && chars[i-1] != ' '), 1:length(chars))

    for i in ind_sign
        sign = chars[i]
        j = i
        while j < length(chars) && (isnumeric(chars[j+1]) || chars[j+1] == ".")
            chars[j] = dict_normal_to_super[chars[j+1]]
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
        rat = stoich_coef_round(num//den)
        replacement = get(dict_frac_unicode, rat, string(rat))
        start_idx = m.offset
        end_idx = start_idx + sizeof(m.match) - 1
        prefix = start_idx > first(eachindex(s)) ? s[first(eachindex(s)):prevind(s, start_idx)] : ""
        suffix = end_idx < sizeof(s) ? s[end_idx+1:end] : ""
        s = prefix * replacement * suffix
    end

    chars = collect(s)

    ind_sign = findall(i -> chars[i] in keys(dict_normal_to_sub) && i > 1 && chars[i-1] != ' ' && !(chars[i-1] in keys(dict_normal_to_sub)), 1:length(chars))

    for i in ind_sign
        j = i
        while j <= length(chars) && chars[j] in keys(dict_normal_to_sub)
            chars[j] = dict_normal_to_sub[chars[j]]
            j += 1
        end
    end

    return join(chars)
end

function unicode_to_phreeqc(s::AbstractString)
    chars = collect(s)

    ind_sign = findall(k -> k == '⁺' || k == '⁻', chars)
    for i in ind_sign   
        sign = chars[i]
        j = i
        while j > 1 && chars[j-1] in keys(dict_super_to_normal)
            chars[j] = chars[j-1]
            j -= 1
        end
        chars[j] = sign
    end
    s = super_to_normal(join(chars))

    s = replace(s, r"[₀₁₂₃₄₅₆₇₈₉]" => c -> dict_sub_to_normal[c[1]])

    pattern = r"(\d+(\.\d+)?|)([¼½¾⅐⅑⅒⅓⅔⅕⅖⅗⅘⅙⅚⅛⅜⅝⅞])"
    matches = collect(eachmatch(pattern, s))
    for m in reverse(matches)
        float_part = isempty(m.captures[1]) ? 0. : parse(Float64, m.captures[1])
        sum_value = float_part + dict_unicode_frac[m.captures[3][1]]
        replacement = string(stoich_coef_round(sum_value))
        start_idx = m.offset
        end_idx = start_idx + sizeof(m.match) - 1
        prefix = start_idx > first(eachindex(s)) ? s[first(eachindex(s)):prevind(s, start_idx)] : ""
        suffix = end_idx < sizeof(s) ? s[end_idx+1:end] : ""
        s = prefix * replacement * suffix
    end

    return s
end

function parse_formula(formula::AbstractString)
    function safe_nextind(s::AbstractString, i::Integer, n::Integer=1)
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

    counts = Dict{Symbol,Number}()

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
            factor = (m === nothing) ? 1 : begin s = m.match ; factor = occursin("//", s) ? parse(Rational{Int}, s) : parse(Float64, s) end
            offset = (m === nothing) ? 0 : length(m.match)

            for (el, n) in parse_formula(inner)
                counts[el] = get(counts, el, 0) + n * factor
            end

            i = safe_nextind(formula, j, offset)
        else
            m = match(r"^(\p{Lu}[\p{Ll}\u0300-\u036F]?)(([0-9]+//[0-9]+)|([0-9]+(?:\.[0-9]+)?))?", formula[i:end])
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

    # T = promote_type(typeof.(stoich_coef_round.(values(counts)))...)

    return Dict(k => stoich_coef_round(v) for (k, v) in counts)

end

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

calculate_molar_mass(atoms::AbstractDict{Symbol,T}) where {T<:Number} =
    uconvert(g/mol, sum(cnt * elements[element].atomic_mass for (element, cnt) in atoms; init=0u) * AvogadroConstant)

function replace_graphemes(s::AbstractString, old_new::Pair...)
    gs = collect(graphemes(s))

    mapping = Dict{String,String}()
    for pair in old_new
        mapping[string(pair.first)] = string(pair.second)
    end

    for i in eachindex(gs)
        if haskey(mapping, gs[i])
            gs[i] = mapping[gs[i]]
        end
    end
    return join(gs)
end

function merge_sum_dicts(dicts::Vector{Dict{Symbol, <:Number}})
    result = Dict{Symbol, Number}()
    for d in dicts
        for (k, v) in d
            result[k] = get(result, k, 0) + v
        end
    end
    return Dict(k => stoich_coef_round(v) for (k, v) in result)
end

function to_mendeleev(oxides::AbstractDict{Symbol,T}) where {T<:Number}
    result = Dict{Symbol, Number}()
    for (ox, coef) in oxides
        if ox ∉ [:Zz, :Zz⁺, :e, :e⁻]
            idx = findfirst(p -> p.first == ox, cement_to_mendeleev)
            idx !== nothing || error("$(ox) is not a valid oxide identifier")
            mend = cement_to_mendeleev[idx].second
            for (k, v) in mend
                result[k] = get(result, k, 0) + v*coef
            end
        end
    end
    return length(result) > 0 ? Dict(k => stoich_coef_round(v) for (k, v) in result) : result
end

const fwd_arrows = ['>', '→', '↣', '↦', '⇾', '⟶', '⟼', '⥟', '⥟', '⇀', '⇁', '⇒', '⟾']
const bwd_arrows = ['<', '←', '↢', '↤', '⇽', '⟵', '⟻', '⥚', '⥞', '↼', '↽', '⇐', '⟽']
const double_arrows = ['↔', '⟷', '⇄', '⇆', '⇌', '⇋', '⇔', '⟺']
const pure_rate_arrows = ['⇐', '⟽', '⇒', '⟾', '⇔', '⟺']
const equal_signs = ['=', '≔', '⩴', '≕']
const EQUAL_REACTION = vcat(fwd_arrows, bwd_arrows, double_arrows, pure_rate_arrows, equal_signs)
const EQUAL_REACTION_SET = Set(EQUAL_REACTION)

function parse_equation(equation::AbstractString)
    equal_sign  = nothing
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
        terms = split(side, " +")  # Split input string by " +" to get each term separately
        result = Dict{String, Real}()  # Initialize dictionary to store formula => coefficient
        
        for term in terms
            t = strip(term)  # Remove leading/trailing whitespace from each term
            
            # Regex to capture optional coefficient (integer, floating or rational a//b) at start, then formula
            m = match(r"^(\d+//\d+|\d*\.?\d+)?\s*(.+)$", t)
            
            if m !== nothing
                coeff_str = m.captures[1]  # The coefficient substring, may be missing
                formula = strip(m.captures[2])    # The chemical formula substring
                
                # If coefficient string is missing or empty, set coefficient to 1
                coeff = coeff_str === nothing || coeff_str == "" ? 1 : eval(Meta.parse(coeff_str))
                
                # Check that coefficient is a real number (Int, Float64, Rational...)
                if !(coeff isa Real)
                    error("Invalid coefficient: $coeff_str")
                end
                
                # Store formula and coefficient in dictionary
                result[formula] = coeff
            else
                error("Unexpected term format: $term")
            end
        end
        
        return Dict(k => stoich_coef_round(v) for (k, v) in result)  # Return dictionary of formula => coefficient
    end
    reactants = left_side == "∅" || left_side == "" ? Dict{String, Int}() : parse_side(left_side)
    products  = right_side == "∅" || right_side == "" ? Dict{String, Int}() : parse_side(right_side)
    return reactants, equal_sign, products
end

function format_equation(coeffs::AbstractDict; scaling=1)
    # Separate reactants and products
    reactants = String[]
    products = String[]
    total_charge_left = 0.0
    total_charge_right = 0.0

    for (species, coeff) in coeffs
        if species !== "Zz"
            coeff = stoich_coef_round(coeff*scaling)

            # Format the coefficient
            abs_coeff = abs(coeff)
            coeff_str = if isapprox(abs_coeff, 1, atol=1e-6)
                ""
            elseif isinteger(abs_coeff)
                string(Int(abs_coeff))
            else
                string(abs_coeff)
            end

            if coeff < 0
                push!(reactants, "$coeff_str$species")
                total_charge_left += coeff * extract_charge(species)
            elseif coeff > 0
                push!(products, "$coeff_str$species")
                total_charge_right += coeff * extract_charge(species)
            end
        end
    end

    # Build the initial equation
    left_side = join(reactants, " + ")
    right_side = join(products, " + ")
    equation = "$left_side = $right_side"

    # Compute the charge difference (corrected)
    charge_diff = total_charge_right + total_charge_left

    # Balance charges if necessary
    if !isapprox(charge_diff, 0, atol=1e-6)
        needed_e = stoich_coef_round(abs(charge_diff))
        e_term = needed_e == 1 ? "e⁻" : "$needed_e" * "e⁻"

        if charge_diff < 0
            # Add e- to the left (reactants)
            left_side = isempty(left_side) ? e_term : "$left_side + $e_term"
        else
            # Add e- to the right (products)
            right_side = isempty(right_side) ? e_term : "$right_side + $e_term"
        end
        equation = "$left_side = $right_side"
    end

    return equation
end

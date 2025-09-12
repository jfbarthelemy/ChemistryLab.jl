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
            factor = (m === nothing) ? 1.0 : begin s = m.match ; factor = occursin("//", s) ? parse(Rational{Int}, s) : parse(Float64, s) end
            offset = (m === nothing) ? 0 : length(m.match)

            for (el, n) in parse_formula(inner)
                counts[el] = get(counts, el, 0) + n * factor
            end

            i = safe_nextind(formula, j, offset)
        else
            m = match(r"^(\p{Lu}\p{Ll}?)([0-9]+//[0-9]+|[0-9]+(?:\.[0-9]+)?)?", formula[i:end])
            if m !== nothing
                el, countstr = m.captures
                el = Symbol(el)
                cnt = if countstr === nothing || isempty(countstr)
                    1
                elseif occursin("//", countstr)
                    parse(Rational{Int}, countstr)
                else
                    parse(Float64, countstr)
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

    T = promote_type(typeof.(stoich_coef_round.(values(counts)))...)

    return Dict(k => convert(T, v) for (k, v) in counts)

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

calculate_molar_mass(atoms::AbstractDict{Symbol,T}) where {T} =
    uconvert(g/mol, sum(cnt * elements[element].atomic_mass for (element, cnt) in atoms) * AvogadroConstant)
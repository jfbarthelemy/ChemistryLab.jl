abstract type Callable end

lib_func(::Val{α}) where {α} = T -> T^α
lib_func(::Val{-1}) = T -> inv(T)
lib_func(::Val{-2}) = T -> inv(T^2)
lib_func(::Val{-3}) = T -> inv(T^3)
lib_func(::Val{0}) = T -> 1
lib_func(::Val{:log}) = T -> log(T/1K)
lib_func(::Val{:logdivT}) = T -> log(T/1K)/T

lib_func_derivative(::Val{α}) where {α} = T -> α * T^(α - 1)
lib_func_derivative(::Val{0}) = T -> 0
lib_func_derivative(::Val{-1}) = T -> -inv(T^2)
lib_func_derivative(::Val{-2}) = T -> -2inv(T^3)
lib_func_derivative(::Val{-3}) = T -> -3inv(T^4)
lib_func_derivative(::Val{:log}) = T -> 1 / T
lib_func_derivative(::Val{:logdivT}) = T -> (1 - log(T/1K)) / T^2

lib_func_primitive(::Val{α}) where {α} = T -> T^(α + 1) / (α + 1)
lib_func_primitive(::Val{-1}) = T -> log(T / 1.0K)
lib_func_primitive(::Val{-2}) = T -> -inv(T)
lib_func_primitive(::Val{-3}) = T -> -inv(T^2)/2
lib_func_primitive(::Val{:log}) = T -> T * (log(T/1K) - 1)
lib_func_primitive(::Val{:logdivT}) = T -> (log(T/1K))^2 / 2

lib_func_double_primitive(::Val{α}) where {α} = T -> T^(α + 2) / ((α+1)*(α+2))
lib_func_double_primitive(::Val{-1}) = T -> (log(T/1K))^2 / 2
lib_func_double_primitive(::Val{-2}) = T -> -log(T/1K)
lib_func_double_primitive(::Val{-3}) = T -> inv(T)/2
lib_func_double_primitive(::Val{:log}) = T -> (T/2)*(log(T/1K))^2 - 2*T*log(T/1K) + (3/2)*T
lib_func_double_primitive(::Val{:logdivT}) = T -> (1/6)*(log(T/1K))^3

lib_func(αs) = (lib_func∘Val).(αs)

const dict_functions = Dict(
    :Cp => (lib_func([0, 1, -2, -0.5, 2, 3, 4, -3, -1, 0.5, :log]), :a, J/(mol*K))
)

function thermo_function(var::Symbol, coeffs::AbstractVector{<: Number}; Tref=298.15K)
    funcs, param, varunit = dict_functions[var]
    var_types = typeof.([varunit/f(1.0K) for f in funcs])
    non_zero_indices = [i-1 for (i, c) in enumerate(coeffs) if i==1 || !iszero(c)]
    type_name = Symbol("Custom$(var)_$(join(string.(non_zero_indices), '_'))")
    fields = [Symbol(param, "$i") for i in non_zero_indices]
    field_decls = [:(($(f)::$(var_types[i+1]))) for (i, f) in zip(non_zero_indices, fields)]

    return @eval begin
        struct $type_name <: Callable
            $(field_decls...)
            Tref::typeof(1.0K)
            function $type_name($(fields...); Tref=$Tref)
                new($(fields...), Tref)
            end
        end

        function (X::$type_name)(T)
            $(Expr(:call, :+,
                [:(X.$(Symbol(param, "$i")) * $funcs[$(i+1)](T))
                 for i in non_zero_indices]...))
        end

        function Base.show(io::IO, X::$type_name)
            # print(io, $(string(var)), "[", $(join(string.(non_zero_indices), ",")), "](")
            print(io, $(string(var)), "(T) with {")
            first = true
            for i in $non_zero_indices
                field = Symbol($(QuoteNode(param)), string(i))
                if !first
                    print(io, ", ")
                end
                print(io, field, "=", getfield(X, field))
                first = false
            end
            print(io, "; Tref=", X.Tref, "}")
        end

        $type_name($(coeffs[non_zero_indices .+ 1]...); Tref=$Tref)
    end
end

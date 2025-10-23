abstract type Callable end  # Base type for all callable thermodynamic functions

# Base thermodynamic functions with optimized implementations for common cases
lib_func(::Val{α}) where {α} = T -> T^α  # General power function
lib_func(::Val{-1}) = T -> inv(T)        # Optimized for T^-1
lib_func(::Val{-2}) = T -> inv(T^2)      # Optimized for T^-2
lib_func(::Val{-3}) = T -> inv(T^3)      # Optimized for T^-3
lib_func(::Val{0}) = T -> 1.0            # Constant function
lib_func(::Val{0.5}) = T -> sqrt(T)      # Optimized square root
lib_func(::Val{-0.5}) = T -> inv(sqrt(T)) # Optimized inverse square root
lib_func(::Val{:log}) = T -> log(ustrip(T))          # Logarithm function
lib_func(::Val{:logdivT}) = T -> log(ustrip(T))/T    # Logarithm divided by T

# Derivative functions
lib_func_derivative(::Val{α}) where {α} = T -> α * T^(α - 1)  # General derivative
lib_func_derivative(::Val{0}) = T -> 0.0                 # Derivative of constant
lib_func_derivative(::Val{-1}) = T -> -inv(T^2)          # Derivative of T^-1
lib_func_derivative(::Val{-2}) = T -> -2inv(T^3)         # Derivative of T^-2
lib_func_derivative(::Val{-3}) = T -> -3inv(T^4)         # Derivative of T^-3
lib_func_derivative(::Val{:log}) = T -> 1 / T            # Derivative of log(T)
lib_func_derivative(::Val{:logdivT}) = T -> (1 - log(ustrip(T))) / T^2  # Derivative of log(T)/T

# Primitive (integral) functions
lib_func_primitive(::Val{α}) where {α} = T -> T^(α + 1) / (α + 1)  # General integral
lib_func_primitive(::Val{-1}) = T -> log(ustrip(T))               # Integral of T^-1
lib_func_primitive(::Val{-2}) = T -> -inv(T)                      # Integral of T^-2
lib_func_primitive(::Val{-3}) = T -> -inv(T^2)/2                  # Integral of T^-3
lib_func_primitive(::Val{:log}) = T -> T * (log(ustrip(T)) - 1)    # Integral of log(T)
lib_func_primitive(::Val{:logdivT}) = T -> (log(ustrip(T)))^2 / 2  # Integral of log(T)/T

# Double primitive (second integral) functions
lib_func_double_primitive(::Val{α}) where {α} = T -> T^(α + 2) / ((α+1)*(α+2))  # General second integral
lib_func_double_primitive(::Val{-1}) = T -> T * (log(ustrip(T)) - 1)            # Second integral of T^-1
lib_func_double_primitive(::Val{-2}) = T -> -log(ustrip(T))                      # Second integral of T^-2
lib_func_double_primitive(::Val{-3}) = T -> inv(T)/2                             # Second integral of T^-3
lib_func_double_primitive(::Val{:log}) = T -> T^2*log(ustrip(T))/2-3T^2/4        # Second integral of log(T)
lib_func_double_primitive(::Val{:logdivT}) = T -> T*(log(ustrip(T)))^2/2-T*log(ustrip(T))+T  # Second integral of log(T)/T

# Pre-compile common function sets for better performance
const Cp_func = (lib_func∘Val).([0, 1, -2, -0.5, 2, 3, 4, -3, -1, 0.5, :log])
const H_func = (lib_func_primitive∘Val).([0, 1, -2, -0.5, 2, 3, 4, -3, -1, 0.5, :log])
const S_func = (lib_func_primitive∘Val).([-1, 0, -3, -1.5, 1, 2, 3, -4, -2, -0.5, :logdivT])
const G_func = (lib_func_double_primitive∘Val).([-1, 0, -3, -1.5, 1, 2, 3, -4, -2, -0.5, :logdivT])

# Dictionary mapping thermodynamic properties to their function sets, parameter names, and units
const dict_functions = Dict(
    :Cp => (Cp_func, :a, J/(mol*K)),  # Heat capacity functions
    :H => (vcat([lib_func(Val(0))], H_func), :aH, J/mol),  # Enthalpy functions (integrals of Cp)
    :S => (vcat([lib_func(Val(0))], S_func), :aS, J/(mol*K)),  # Entropy functions
    :G => (vcat((lib_func∘Val).([0, 1]), G_func), :aG, J/mol),  # Gibbs energy functions (double integrals)
)

# Cache for generated types to avoid redefinition
const ThermoTypeCache = Dict{Tuple, DataType}()

"""
    generate_thermo_type(var::Symbol, non_zero_indices::Vector{Int}, with_units::Bool, Tref)

Generate a unique callable struct type for the thermodynamic property `var` with given non-zero coefficient indices.

This function performs dynamic type and method generation at global module level,
and caches the generated type for future construction.
"""
function generate_thermo_type(var::Symbol, non_zero_indices::Vector{Int}, default_types::Vector{DataType}; with_units=true, Tref=with_units ? 298.15K : 298.15, startindex=0)
    key = (var, tuple(non_zero_indices...), with_units, startindex)
    if haskey(ThermoTypeCache, key)
        return ThermoTypeCache[key]
    end

    funcs, param, varunit = dict_functions[var]
    var_types = with_units ? typeof.([varunit / f(1.0K) for f in funcs]) : default_types

    fields = [Symbol(param, "$i") for i in non_zero_indices]
    field_decls = [:(($(f)::$(var_types[i+1-startindex]))) for (i, f) in zip(non_zero_indices, fields)]

    type_name = Symbol("Custom$(var)_$(join(string.(non_zero_indices), '_'))" * (with_units ? "" : "_nounit"))

    # Define struct
    type_expr = :(
        struct $type_name <: Callable
            $(field_decls...)
            Tref::typeof($Tref)
            function $type_name($(fields...); Tref=$Tref)
                new($(fields...), Tref)
            end
        end
    )

    # Define call method optimized for single or multiple coefficients
    if length(non_zero_indices) == 1
        call_expr = :(
            function (X::$type_name)(T)
                @inbounds X.$(Symbol(param, "$(non_zero_indices[1])")) * $(funcs[non_zero_indices[1]+1-startindex])(T)
            end
        )
    else
        call_expr = :(
            function (X::$type_name)(T)
                    @inbounds $(Expr(:call, :+,
                        [:(X.$(Symbol(param, "$i")) * $(funcs[i+1-startindex])(T))
                         for i in non_zero_indices]...))
            end
        )
    end

    # Define show method
    show_exprs = Expr[]
    for (idx, i) in enumerate(non_zero_indices)
        field = Symbol(param, "$i")
        push!(show_exprs, :(print(io, $(QuoteNode(field)), "=", getfield(X, $(QuoteNode(field))))))
        if idx < length(non_zero_indices)
            push!(show_exprs, :(print(io, ", ")))
        end
    end

    show_expr = :(
        function Base.show(io::IO, X::$type_name)
            print(io, $(string(var)), "(T) with {")
            $(show_exprs...)
            print(io, "; Tref=", X.Tref, "}")
        end
    )

    # Define method call without argument (uses reference temperature)
    call_no_arg_expr = :(
        function (X::$type_name)()
            X(X.Tref)
        end
    )

    # Evaluate expressions *globally* in the module to avoid world-age issues
    Base.eval(@__MODULE__, type_expr)
    Base.eval(@__MODULE__, call_expr)
    Base.eval(@__MODULE__, show_expr)
    Base.eval(@__MODULE__, call_no_arg_expr)

    # Cache and return the generated type, using invokelatest to avoid world age issues when accessing binding
    ThermoTypeCache[key] = Base.invokelatest(getfield, @__MODULE__, type_name)
    return ThermoTypeCache[key]
end

"""
    thermo_function(var::Symbol, coeffs::AbstractVector{<:Number}; with_units=true, Tref=298.15K, startindex=0)

Create an instance of callable thermodynamic function for property `var` with specified coefficients.

This function calls `generate_thermo_type` once per unique coefficient pattern,
then quickly constructs and returns an instance of the generated type.
"""
function thermo_function(var::Symbol, coeffs::AbstractVector{<:Number}; with_units=true, Tref=with_units ? 298.15K : 298.15, startindex=0)
    # Identify indices of non-zero coefficients (always include the first)
    non_zero_indices = [i-1+startindex for (i,c) in enumerate(coeffs) if i == 1 || !iszero(c)]
    non_zero_coeffs = [coeffs[i+1-startindex] for i in non_zero_indices]

    # Generate or reuse the callable struct type for these coefficients
    type_t = generate_thermo_type(var, non_zero_indices, typeof.(coeffs); with_units=with_units, Tref=Tref, startindex=startindex)

    # Construct and return an instance with given coefficients and reference temperature
    return Base.invokelatest(type_t, non_zero_coeffs...; Tref=Tref)
end


# Example usage:
# cp = thermo_function(:Cp, [1.0, 0.0, 2.0, 0.0, 3.0])
# val_at_300 = cp(300.0)
# println(cp)

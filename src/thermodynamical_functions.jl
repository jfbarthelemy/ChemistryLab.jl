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

struct ThermoFunction{F<:NamedTuple,C<:NamedTuple,T<:Number,Z<:Number} <: Callable
    name::Symbol
    bases::F
    coeffs::C
    Tref::T
    zeroinit::Z
end

function (lf::ThermoFunction)(T)
    s = lf.zeroinit
    @inbounds @simd for name in keys(lf.coeffs)
        s += getfield(lf.coeffs, name) * getfield(lf.bases, name)(T)
    end
    return s
end

function (lf::ThermoFunction)()
    return lf(lf.Tref)
end

coefficients(lf::ThermoFunction) = lf.coeffs

function ThermoFunction(name::Symbol, coeffs::AbstractVector{<:Number}; Tref=298.15, startindex=0)

    with_units = promote_type(typeof.(coeffs)...) <: Quantity
    if with_units && !(Tref isa Quantity)
        Tref *= K
    end

    funcs, param, varunit = dict_functions[name]

    nonzero = findall(!iszero, coeffs)
    kept_names = Symbol.(string.(param), string.("_", nonzero .+ (startindex - 1)))
    kept_funcs = funcs[nonzero]
    kept_coefs = coeffs[nonzero]

    base_nt = NamedTuple{Tuple(kept_names)}(Tuple(kept_funcs))
    coef_nt = NamedTuple{Tuple(kept_names)}(Tuple(kept_coefs))

    return ThermoFunction(name, base_nt, coef_nt, Tref, with_units ? 0*varunit : 0)
end

function Base.show(io::IO, lf::ThermoFunction)
    print(io, lf.name, "(T) with ")
    print(io, lf.coeffs)
    print(io, " and Tref = ", lf.Tref)
end

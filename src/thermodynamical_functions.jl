abstract type Callable end

struct Cp <: Callable
    a0::typeof(1.0J/K/mol)
    a1::typeof(1.0J/K^2/mol)
    a2::typeof(1.0J*K/mol)
    a3::typeof(1.0J/√K/mol)
end
(cp::Cp)(T) = cp.a0 + cp.a1 * T + cp.a2 / T^2 + cp.a3 / √T

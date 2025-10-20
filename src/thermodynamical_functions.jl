abstract type Callable end

const Cp_types = [
        typeof(1.0*J/(mol*K)),
        typeof(1.0*J/(mol*K^2)),
        typeof(1.0*(J*K)/mol),
        typeof(1.0*J/(mol*K^0.5)),
        typeof(1.0*J/(mol*K^3)),
        typeof(1.0*J/(mol*K^4)),
        typeof(1.0*J/(mol*K^5)),
        typeof(1.0*(J*K^2)/mol),
        typeof(1.0*J/mol),
        typeof(1.0*J/(mol*K^1.5)),
        typeof(1.0*J/(mol*K)),
    ]

struct Cp <: Callable
    a0::typeof(1.0*J/(mol*K))
    a1::typeof(1.0*J/(mol*K^2))
    a2::typeof(1.0*(J*K)/mol)
    a3::typeof(1.0*J/(mol*K^0.5))
    a4::typeof(1.0*J/(mol*K^3))
    a5::typeof(1.0*J/(mol*K^4))
    a6::typeof(1.0*J/(mol*K^5))
    a7::typeof(1.0*(J*K^2)/mol)
    a8::typeof(1.0*J/mol)
    a9::typeof(1.0*J/(mol*K^1.5))
    a10::typeof(1.0*J/(mol*K))
end

function Cp(args...)
    n = length(Cp_types)
    vals = ntuple(i -> i <= length(args) ? args[i] : zero(Cp_types[i]), n)
    return Cp(vals...)
end

(cp::Cp)(T) = cp.a0 + cp.a1*T + cp.a2/T^2 + cp.a3/√T + cp.a4*T^2 + cp.a5*T^3 + cp.a6*T^4 + cp.a7/T^3 + cp.a8/T + cp.a9*√T + cp.a10*log(T/1.0K)

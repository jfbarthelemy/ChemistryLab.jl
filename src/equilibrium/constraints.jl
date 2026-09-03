# ── constraints.jl ────────────────────────────────────────────────────────────
#
# Equilibrium under something other than fixed (T, P).
#
# The vehicle is the solver's parameter block: a prescribed property adds one
# unknown — the temperature, the pressure — and one equation, its own residual, to
# the SAME square system that carries the amounts and the multipliers. There is no
# outer loop around the equilibrium solve, which is the structure Reaktoro uses
# and the reason an adiabatic solve costs one extra equation rather than a second
# solve per trial temperature.
#
# Reaktoro has a second vehicle, the implicit titrant, for a prescribed chemical
# potential — pH, pE, a fixed fugacity. That one adds a COLUMN to the conservation
# matrix rather than a parameter, and is not implemented here. Asking for it
# raises rather than silently ignoring it.

"""
    EquilibriumConstraint

What is held fixed while the Gibbs energy is minimized. One of
[`FixedTP`](@ref), [`FixedEnthalpy`](@ref), [`Adiabatic`](@ref),
[`FixedVolume`](@ref) or [`SealedVolume`](@ref).
"""
abstract type EquilibriumConstraint end

"""
    FixedTP()

Temperature and pressure both held at the values the state carries. The default,
and the only constraint the interior-point back ends can honor.
"""
struct FixedTP <: EquilibriumConstraint end

"""
    FixedEnthalpy(H)

Enthalpy prescribed, temperature unknown. `H` is a quantity with the dimensions
of energy, referred to the same element basis as [`enthalpy`](@ref) — so in
practice it comes from `enthalpy(state)` of another state on the same budget,
never from a table.
"""
struct FixedEnthalpy{Q} <: EquilibriumConstraint
    H::Q
end

"""
    Adiabatic()

Enthalpy conserved at the value the initial state carries, temperature unknown.
This is a reaction in a vessel that exchanges no heat: the enthalpy released by
the reaction goes into raising the temperature.

Equivalent to `FixedEnthalpy(enthalpy(state))`, and preferable because it cannot
be given an enthalpy from a different element basis by mistake.
"""
struct Adiabatic <: EquilibriumConstraint end

"""
    FixedVolume(V)

Volume prescribed, pressure unknown. `V` is a quantity with the dimensions of
volume, and only the species carrying a molar volume contribute — see
[`volume`](@ref).

!!! warning "Needs a system whose volume depends on pressure"
    A condensed system's does not. The molar volumes of water and of the minerals
    in the shipped databases are exactly pressure-independent, so prescribing the
    volume of a paste or of an aqueous solution is an equation the pressure cannot
    satisfy, and the constructor refuses it rather than diverging — the error
    names the lever it measured. Declare a gas phase to make the constraint
    meaningful. To report the volume change of a sealed specimen at fixed
    pressure, which is what a hydrating binder needs, use
    `porosity(state, reference)` instead.
"""
struct FixedVolume{Q} <: EquilibriumConstraint
    V::Q
end

"""
    SealedVolume()

Volume held at the value the initial state occupies, pressure unknown: a rigid
sealed vessel. Equivalent to `FixedVolume(volume(state).total)`.
"""
struct SealedVolume <: EquilibriumConstraint end

"""
    _total_enthalpy(system, n, T, P)

`Σᵢ nᵢ ΔₐH⁰ᵢ(T, P)` in joules, as a bare number, over the species that carry an
enthalpy of formation. The element type follows `n`, `T` and `P`, so this
differentiates.
"""
function _total_enthalpy(system::ChemicalSystem, n, T, P)
    tot = zero(promote_type(eltype(n), typeof(T), typeof(P)))
    for (i, s) in enumerate(system.species)
        f = _std_property(s, :ΔₐH⁰)
        f === nothing && continue
        tot += n[i] * ustrip(us"J/mol", f(T = T * u"K", P = P * u"Pa"; unit = true))
    end
    return tot
end

"""
    _total_volume(system, n, T, P)

`Σᵢ nᵢ V⁰ᵢ(T, P)` in cubic meters, as a bare number, over the species that carry
a molar volume.
"""
function _total_volume(system::ChemicalSystem, n, T, P)
    tot = zero(promote_type(eltype(n), typeof(T), typeof(P)))
    for (i, s) in enumerate(system.species)
        _has_molar_volume(s) || continue
        tot += n[i] * ustrip(
            us"m^3/mol", _molar_volume(s)(T = T * u"K", P = P * u"Pa"; unit = true),
        )
    end
    return tot
end

"""
    _constraint_blocks(constraint, des, state, p, n0) -> NamedTuple

The parameter block a constraint contributes: `nq`, the callbacks `gq`, `hq` and
`cq`, the starting guess `q0`, the difference-step scale `qscale`, and `apply`,
which writes the parameter that was found back onto the resulting state.

`nq = 0` for [`FixedTP`](@ref), in which case the solver takes its plain path.
"""
function _constraint_blocks(::FixedTP, des, state, p, n0)
    return (;
        nq = 0, gq = nothing, hq = nothing, cq = nothing,
        Aq = zeros(Float64, size(des === nothing ? zeros(0, 0) : des.A, 1), 0),
        q0 = Float64[], qscale = Float64[],
        apply = (T, P, q) -> (T, P),
    )
end

# Temperature unknown: `q[1]` is T in kelvin.
function _temperature_blocks(des, state, p, target_H)
    system = state.system
    P = p.P
    # The reference potentials and the activity model are both functions of T, so
    # both are rebuilt at the current guess. Leaving the activity model at the
    # starting temperature would minimize a Gibbs energy that does not exist.
    gq = (q, params) -> [
        ustrip(
            us"J/mol",
            s[:ΔₐG⁰](T = q[1] * u"K", P = P * u"Pa"; unit = true),
        ) / (ustrip(us"J/(mol*K)", Constants.R) * q[1])
            for s in system.species
    ]
    hq = (x, q, params) -> des.lna(x, merge(params, (T = q[1],)))
    # The residual is scaled by RT so it is dimensionless and comparable with the
    # stationarity rows, which are in RT units. Unscaled it is 10⁵ J and swamps
    # every other row of the Newton system.
    cq = (x, q, params) -> [
        (_total_enthalpy(system, x, q[1], P) - target_H) /
            (ustrip(us"J/(mol*K)", Constants.R) * q[1] * max(sum(x), 1.0)),
    ]
    return (;
        nq = 1, gq = gq, hq = hq, cq = cq,
        Aq = zeros(Float64, size(des.A, 1), 1),
        q0 = [p.T], qscale = [p.T],
        apply = (T, Pv, q) -> (q[1] * u"K", Pv),
    )
end

"""
    _pressure_lever(system, n, T, P) -> Real

Relative sensitivity of the system's volume to pressure, `(∂V/∂P)·P/V`, by a
central difference over one percent of `P`.

A volume constraint prescribes `V` and solves for `P`, so it needs this to be
non-negligible. It usually is not for a condensed system: in the databases
shipped here the molar volumes of water and of the minerals do not depend on
pressure at all — `V⁰(1 bar) = V⁰(100 bar)` exactly for `H2O@` and `Cal` — and
only the partial molar volumes of a few aqueous ions vary, `OH-` by 8 % over
100 bar. A kilogram of water with 0.25 mol of ions then has a lever of about
`1e-6`, meaning some 10 000 bar to change the volume by one percent. Newton on
that residual takes an enormous step and the pressure leaves the domain of the
equation of state.

Which is the physics, not a solver defect: the volume of an incompressible
condensed system is fixed by its composition, and pressure has no purchase on it.
A gas phase gives it one.
"""
function _pressure_lever(system::ChemicalSystem, n, T, P)
    δ = 0.01 * P
    Vp = _total_volume(system, n, T, P + δ)
    Vm = _total_volume(system, n, T, P - δ)
    V = _total_volume(system, n, T, P)
    V == 0 && return 0.0
    return abs((Vp - Vm) / (2δ)) * P / abs(V)
end

# Pressure unknown: `q[1]` is P in pascals.
function _pressure_blocks(des, state, p, target_V)
    system = state.system
    T = p.T

    n_now = Float64[ustrip(us"mol", x) for x in state.n]
    lever = _pressure_lever(system, n_now, T, p.P)
    if lever < 1.0e-4
        throw(
            ArgumentError(
                "a volume constraint cannot be solved on this system: its volume " *
                    "barely depends on pressure, the relative lever (∂V/∂P)·P/V " *
                    "being $(round(lever, sigdigits = 2)). Changing the volume by " *
                    "one percent would take about " *
                    "$(round(0.01 / max(lever, 1.0e-30) * p.P / 1.0e5, sigdigits = 2)) " *
                    "bar. The molar volumes of water and of the minerals in the " *
                    "shipped databases do not depend on pressure at all, so a " *
                    "condensed system's volume is fixed by its composition and " *
                    "there is nothing for the pressure to do. Declare a gas phase, " *
                    "or hold the pressure and let the volume follow — " *
                    "`porosity(state, reference)` reports the volume change under " *
                    "sealed curing without needing a constraint.",
            )
        )
    end
    gq = (q, params) -> [
        ustrip(
            us"J/mol",
            s[:ΔₐG⁰](T = T * u"K", P = q[1] * u"Pa"; unit = true),
        ) / (ustrip(us"J/(mol*K)", Constants.R) * T)
            for s in system.species
    ]
    hq = (x, q, params) -> des.lna(x, merge(params, (P = q[1],)))
    # Scaled by the target volume: the residual is then a relative volume error.
    cq = (x, q, params) -> [
        (_total_volume(system, x, T, q[1]) - target_V) / max(abs(target_V), 1.0e-12),
    ]
    return (;
        nq = 1, gq = gq, hq = hq, cq = cq,
        Aq = zeros(Float64, size(des.A, 1), 1),
        q0 = [p.P], qscale = [max(p.P, 1.0e5)],
        apply = (Tv, P, q) -> (Tv, q[1] * u"Pa"),
    )
end

_constraint_blocks(c::FixedEnthalpy, des, state, p, n0) =
    _temperature_blocks(des, state, p, ustrip(us"J", c.H))

_constraint_blocks(::Adiabatic, des, state, p, n0) =
    _temperature_blocks(des, state, p, _total_enthalpy(state.system, n0, p.T, p.P))

_constraint_blocks(c::FixedVolume, des, state, p, n0) =
    _pressure_blocks(des, state, p, ustrip(us"m^3", c.V))

_constraint_blocks(::SealedVolume, des, state, p, n0) =
    _pressure_blocks(des, state, p, _total_volume(state.system, n0, p.T, p.P))

# ── the implicit titrant: a prescribed chemical potential ────────────────────
#
# The second vehicle. A prescribed activity or pH does not add a parameter: it
# adds a COLUMN to the conservation matrix — an unknown amount of a titrant
# substance the system may draw on — and one equation, the prescribed potential.
# The system is open to that substance and closed to everything else.
#
# `Aq` carries the column, so the linear rows read `A x − A[:, titrant] q = b`:
# `q` moles of titrant enter the budget. The residual is the prescribed
# `ln aᵢ`, read from the activity model at the current composition.

"""
    FixedActivity(species, a; titrant = species)

Activity of `species` prescribed, the system open to `titrant`.

`titrant` is the substance whose amount adjusts to hold the activity — by default
the species itself. Its amount is an unknown of the same system, so the answer is
both the equilibrium composition and how much titrant it took to get there, which
is what a titration measures.

The element balance is then `A x − A[:, titrant] q = b`: the system is open to
that one substance and closed to every other.
"""
struct FixedActivity{S, R} <: EquilibriumConstraint
    species::S
    a::R
    titrant::S
end
FixedActivity(species, a; titrant = species) = FixedActivity(species, a, titrant)

"""
    FixedpH(pH; titrant = "H+")

pH prescribed, the system open to `titrant`. Shorthand for
`FixedActivity("H+", 10^-pH; titrant = titrant)`.

`titrant = "H+"` means acid is added or removed. To titrate with a base instead,
name it: `FixedpH(12.5; titrant = "OH-")`. The answer carries the amount, so the
titrant is not a numerical device — it is the reagent consumed.

!!! note "Which pH is prescribed"
    This constrains `−log₁₀ a(H⁺)` in the activity model's own convention, and it
    holds it exactly: the residual is that quantity, and it comes out at the
    target to five decimals. [`pH`](@ref) reports something slightly different —
    `−log₁₀ c(H⁺)` with the concentration taken over the computed liquid volume —
    and with [`DiluteSolutionModel`](@ref) the two differ by about **0.0013**,
    which is `log₁₀` of the `ρ ≈ 1 kg/L` approximation that model makes when it
    converts molality to concentration. The gap grows with ionic strength; use
    [`HKFActivityModel`](@ref) or [`DaviesActivityModel`](@ref) where it matters.
"""
struct FixedpH{R, S} <: EquilibriumConstraint
    pH::R
    titrant::S
end
FixedpH(pH; titrant = "H+") = FixedpH(pH, titrant)

"""
    _titrant_blocks(des, state, p, i_species, i_titrant, ln_a_target)

The parameter block of a prescribed chemical potential: one unknown, the titrant
amount; one linear column, its formula; one equation, `ln aᵢ = ln a_target`.
"""
function _titrant_blocks(des, state, p, i_species, i_titrant, ln_a_target)
    # `−A[:, titrant]`: the titrant ADDS to the budget, and the rows are written
    # as `A x + Aq q − b = 0`.
    Aq = reshape(-des.A[:, i_titrant], size(des.A, 1), 1)
    # A mole scale for the difference step, taken from the system's own size so
    # it is neither absurdly large nor below the resolution of the balance.
    scale = max(sum(Float64[ustrip(us"mol", x) for x in state.n]), 1.0) * 1.0e-6
    cq = (x, q, params) -> [des.lna(x, params)[i_species] - ln_a_target]
    return (;
        nq = 1, gq = (q, params) -> params.ΔₐG⁰overT, hq = nothing, cq = cq,
        Aq = Aq, q0 = [0.0], qscale = [scale],
        apply = (T, P, q) -> (T, P),
        titrant_amount = q -> q[1],
    )
end

function _constraint_blocks(c::FixedActivity, des, state, p, n0)
    i_s = _species_index(des, c.species)
    i_t = _species_index(des, c.titrant)
    return _titrant_blocks(des, state, p, i_s, i_t, log(c.a))
end

function _constraint_blocks(c::FixedpH, des, state, p, n0)
    i_s = _species_index(des, "H+")
    i_t = _species_index(des, c.titrant)
    return _titrant_blocks(des, state, p, i_s, i_t, -c.pH * log(10))
end

"""
    _species_index(des, s) -> Int

Position of a species in the solver's system, by symbol or by object.
"""
function _species_index(des, sym::AbstractString)
    i = findfirst(==(sym), symbol.(des.system.species))
    i === nothing && throw(
        ArgumentError(
            "`$sym` is not a species of this system, so it cannot be a " *
                "titrant or carry a prescribed activity."
        )
    )
    return i
end
_species_index(des, s::AbstractSpecies) = _species_index(des, symbol(s))

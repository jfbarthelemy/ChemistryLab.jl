# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

"""
    KineticsOrdinaryDiffEqExt

Extension activated when `OrdinaryDiffEq` is loaded. Provides the concrete
`ChemistryLab.integrate(::KineticsProblem, ::KineticsSolver; ...)` implementation
using `OrdinaryDiffEq.ODEProblem` and registers `Rodas5P()` as the default solver.

Follows the Leal et al. (2017) formulation with partial equilibrium and
optional thermal coupling (isothermal or semi-adiabatic calorimetry).

# Usage

```julia
using ChemistryLab, OrdinaryDiffEq
sol = integrate(kp, KineticsSolver(; ode_solver=Rodas5P()))
# or shortcut (uses default Rodas5P):
sol = integrate(kp)
```
"""
module KineticsOrdinaryDiffEqExt

using OrdinaryDiffEq
using SciMLBase
import ChemistryLab
import ChemistryLab:
    integrate,
    KineticsProblem,
    KineticsSolver,
    build_kinetics_ode,
    build_u0,
    build_kinetics_params,
    _DEFAULT_KINETICS_SOLVER_FACTORY,
    SemiAdiabaticCalorimeter,
    IsothermalCalorimeter,
    respeciate!,
    _with_equilibrium_solver,
    symbol

# ── Concrete integrate implementation ────────────────────────────────────────

"""
    integrate(kp::KineticsProblem, ks::KineticsSolver) -> ODESolution

Integrate the kinetics ODE using `OrdinaryDiffEq` (Leal et al. 2017 formulation).

The ODE function, initial state, and parameters are built from `kp`.
Calorimetry (isothermal or semi-adiabatic) is integrated directly in the ODE
right-hand-side — no separate `extend_ode!` step.

Default tolerances: `reltol = 1e-8`, `abstol = 1e-10`.

# Examples

```julia
using ChemistryLab, OrdinaryDiffEq

ks  = KineticsSolver(; ode_solver=Rodas5P(), reltol=1e-8, abstol=1e-10)
sol = integrate(kp, ks)
```
"""
function integrate(kp::KineticsProblem, ks::KineticsSolver; kwargs...)
    # `KineticsSolver` also carries an equilibrium solver, and its docstring
    # advertises passing one there. Honor it: without this the field is dead
    # and re-speciation silently never happens.
    kp = _with_equilibrium_solver(kp, ks.equilibrium_solver)

    f! = build_kinetics_ode(kp)
    u0 = build_u0(kp)
    p = build_kinetics_params(kp)

    # A semi-adiabatic cell driven by an equilibrium-coupled model would be
    # driven by the wrong heat. The temperature equation takes its source from
    # `heat_rate`, i.e. from the KINETIC reactions, and under partial equilibrium
    # those only dissolve the anhydrous phases into ions — the hydrates are
    # precipitated by the Gibbs minimization, whose heat that sum cannot see.
    # On an ordinary Portland cement this put the temperature rise at 207 K
    # against the few tens of kelvin a Langavant test gives, and nothing in the
    # run said so. Until the source accounts for the precipitation, say it here.
    if kp.calorimeter isa IsothermalCalorimeter && !isnothing(kp.equilibrium_solver)
        @warn """isothermal calorimetry is coupled to an equilibrium solver: the Q state \
        this integrates accumulates the heat of the KINETIC reactions alone, which under \
        partial equilibrium is the heat of DISSOLUTION and omits the precipitation of the \
        hydrates. `cumulative_heat` and `heat_flow` will report that partial figure. Use \
        `heat_release`, which differences the enthalpy of certified speciations."""
    end

    if kp.calorimeter isa SemiAdiabaticCalorimeter && !isnothing(kp.equilibrium_solver)
        @warn """semi-adiabatic calorimetry is coupled to an equilibrium solver: the \
        temperature is driven by the heat of the KINETIC reactions alone, which under \
        partial equilibrium is the heat of DISSOLUTION and omits the precipitation of \
        the hydrates. The temperature will be badly overestimated. Use \
        `IsothermalCalorimeter` with `heat_release`, which reads certified \
        speciations, until the source term accounts for the equilibrium partition."""
    end

    # Warn for missing Cp° when semi-adiabatic
    if kp.calorimeter isa SemiAdiabaticCalorimeter
        missing_cp = String[]
        for (sp, cp_fn) in zip(kp.system.species, p.cp_fns)
            isnothing(cp_fn) && push!(missing_cp, string(symbol(sp)))
        end
        if !isempty(missing_cp)
            shown = join(missing_cp[1:min(5, length(missing_cp))], ", ")
            suffix = length(missing_cp) > 5 ? "…" : ""
            @warn "SemiAdiabaticCalorimeter: variable Cp_total requires Cp° data per " *
                "species. Missing for $(length(missing_cp)) species " *
                "($shown$suffix). Their contribution to Cp_total is treated as zero."
        end
    end

    # Precedence: call-site kwargs beat the solver's, which beat the defaults.
    # `integrate(kp; reltol = …)` forwards through the one-argument method and
    # used to hit a MethodError, this one taking no kwargs at all.
    defaults = (reltol = 1.0e-8, abstol = 1.0e-10)
    merged = merge(merge(defaults, ks.kwargs), values(kwargs))
    solver = isnothing(ks.ode_solver) ? Rodas5P() : ks.ode_solver

    prob = ODEProblem(f!, u0, kp.tspan, p)

    # Operator splitting: the ODE advances the kinetic minerals with the
    # speciation frozen, and this callback re-equilibrates once per accepted
    # step. Nothing in `u` is touched, so `save_positions = (false, false)`.
    if p.n_be > 0
        # A non-converged solve is a warning, not an exception, and its result is
        # used anyway — so it never reached `eq_failures`. Count it over the run.
        nonconv0 = ChemistryLab.NONCONVERGED[]
        respeciate!(p, u0)          # start from an equilibrated state
        cb = DiscreteCallback(
            (u, t, integrator) -> true,
            integrator -> begin
                # Flag the trajectory: only speciations computed here belong to
                # the solution. Everything else is a probe.
                integrator.p.on_accepted[] = true
                respeciate!(integrator.p, integrator.u)
                integrator.p.on_accepted[] = false
                SciMLBase.u_modified!(integrator, false)
            end;
            save_positions = (false, false),
        )
        sol = solve(prob, solver; callback = cb, merged...)

        if p.eq_failures[] > 0
            @warn "re-speciation failed on $(p.eq_failures[]) step(s); those steps kept a frozen composition."
        end
        nonconv = ChemistryLab.NONCONVERGED[] - nonconv0
        if nonconv > 0
            res = p.eq_worst_residual[]
            abs_res = p.eq_worst_abs_acc[]
            abs_all = p.eq_worst_abs[]
            @warn """$nonconv equilibrium solve(s) stopped short of the optimizer's \
            tolerance and were used anyway. Judge them on the element balance, not \
            on the retcode: worst |Aₑn − bₑ|∞ over the run was \
            $(round(abs_res; sigdigits = 3)) mol ON THE ACCEPTED STEPS, i.e. on the \
            trajectory itself; $(round(abs_all; sigdigits = 3)) mol counting also the \
            Jacobian probes and rejected steps, which never enter the solution. \
            Read the first figure: 1e-10 mol is machine precision whatever the \
            system, while 1e-2 mol against a 0.3 mol sulfate budget is not. \
            How much it matters depends on the RATE LAWS: `bₑ` is integrated from \
            the rates alone, so a law that reads only its own degree of reaction \
            (Parrot-Killoh, Waller) gives a trajectory independent of the \
            speciation, and this figure then bears on the reported composition \
            only — recover that with `speciated_states`, which certifies each instant against the KKT conditions. A \
            law reading log-activities (a saturation ratio) does feed the \
            speciation back into the trajectory, and there this figure is a \
            direct measure of the error. \
            Do NOT simply loosen the optimizer tolerance — on the calcite reference \
            case that degrades the speciation from 4 % to 250 % against Reaktoro."""
        end
    else
        sol = solve(prob, solver; merged...)
    end
    _warn_if_unphysical(sol, p, kp)
    return sol
end

"""
    _warn_if_unphysical(sol, p, kp)

Warn when the trajectory ends on amounts no chemistry can produce.

This is not tidying. Measured on calcite dissolving under `r = k(1 − Ω)` over
`10⁵ s`, `Rodas5P` — and `OrdinaryDiffEq`'s default polyalgorithm, which picks a
stiff method here — return a reaction extent of **−457 mol** and
`retcode = Success`, while the explicit `Tsit5` gets the right answer
(`1.104e-4`) in 85 626 steps. The stiff methods need a Jacobian, the residual
carries a re-speciation the Jacobian does not see, and the error control built on
that same Jacobian reports success on nonsense.

The implicit step of `kinetic_step` does not have the problem — it is exact on the
same case at `Δt = 10³ s` and, at `10⁵ s`, wrong but **reporting** it through its
certificate — and `kinetic_step_adaptive` reaches the equilibrium values to eight
digits in seven steps. Until the ODE route carries the exact Jacobian, an answer
from it deserves this check.
"""
function _warn_if_unphysical(sol, p, kp)
    u = sol.u[end]
    nb, nk, nr = p.n_be, p.n_nk, p.n_rxn_state
    # The test is CREATION OF MATTER, not the sign of an extent. A negative
    # extent is legitimate — that is precipitation — and on the measured failure
    # the extent came out at −457 mol, so a threshold on it would have to be
    # arbitrary. What cannot happen is an amount exceeding what the system was
    # given: `Rodas5P` returned 457 mol of calcite from a budget of 55.6.
    total0 = sum(p.n_initial_full)
    ceiling = 2 * total0 + 1
    bad = String[]
    for j in 1:nk
        v = u[nb + j]
        name = symbol(kp.system.species[kp.idx_kinetic[j]])
        v < -1.0e-8 && push!(bad, "$name = $(round(v, sigdigits = 4)) mol, negative")
        v > ceiling && push!(
            bad,
            "$name = $(round(v, sigdigits = 4)) mol against a total budget of " *
                "$(round(total0, sigdigits = 4))",
        )
    end
    for j in 1:nr
        ξ = u[nb + nk + j]
        abs(ξ) > ceiling && push!(bad, "extent $j = $(round(ξ, sigdigits = 4)) mol")
    end
    isempty(bad) && return nothing
    @warn """the trajectory ends on amounts no chemistry can produce, and the     integrator reported success: $(join(bad, ", ")). A stiff method here needs a     Jacobian that the re-speciation in the residual does not supply, so its error     control is built on the wrong derivative. Use `kinetic_step_adaptive`, whose     step is chosen from a local error estimate on the extents and which is exact on     this class of problem, or an explicit method if the system allows it.""" maxlog = 1
    return nothing
end

# ── __init__: register default solver ────────────────────────────────────────

function __init__()
    return _DEFAULT_KINETICS_SOLVER_FACTORY[] =
        () -> KineticsSolver(; ode_solver = Rodas5P())
end

end  # module KineticsOrdinaryDiffEqExt

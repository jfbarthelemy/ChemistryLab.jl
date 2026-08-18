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
    return sol
end

# ── __init__: register default solver ────────────────────────────────────────

function __init__()
    return _DEFAULT_KINETICS_SOLVER_FACTORY[] =
        () -> KineticsSolver(; ode_solver = Rodas5P())
end

end  # module KineticsOrdinaryDiffEqExt

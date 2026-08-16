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
function integrate(kp::KineticsProblem, ks::KineticsSolver)
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

    defaults = (reltol = 1.0e-8, abstol = 1.0e-10)
    merged = merge(defaults, ks.kwargs)
    solver = isnothing(ks.ode_solver) ? Rodas5P() : ks.ode_solver

    prob = ODEProblem(f!, u0, kp.tspan, p)

    # Operator splitting: the ODE advances the kinetic minerals with the
    # speciation frozen, and this callback re-equilibrates once per accepted
    # step. Nothing in `u` is touched, so `save_positions = (false, false)`.
    if p.n_be > 0
        respeciate!(p, u0)          # start from an equilibrated state
        cb = DiscreteCallback(
            (u, t, integrator) -> true,
            integrator -> (respeciate!(integrator.p, integrator.u); SciMLBase.u_modified!(integrator, false));
            save_positions = (false, false),
        )
        sol = solve(prob, solver; callback = cb, merged...)
    else
        sol = solve(prob, solver; merged...)
    end

    if p.n_be > 0 && p.eq_failures[] > 0
        @warn "re-speciation failed on $(p.eq_failures[]) step(s); those steps kept a frozen composition."
    end
    return sol
end

# ── __init__: register default solver ────────────────────────────────────────

function __init__()
    return _DEFAULT_KINETICS_SOLVER_FACTORY[] =
        () -> KineticsSolver(; ode_solver = Rodas5P())
end

end  # module KineticsOrdinaryDiffEqExt

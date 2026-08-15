# Oracle for the kinetics/equilibrium coupling (test/coupling_reference.jl).
#
# Calcite dissolves at a *constant* rate, so the kinetic half of Leal's system
# is analytic:
#
#     n_Cal(t) = n_Cal(0) - k t,      b_e(t) = b_e(0) + k t * (one CaCO3)
#
# and the whole content of the coupling reduces to whether the equilibrium
# partition follows n_e(t) = phi(b_e(t)), Leal et al. (2017) Eq. 54. Fixing the
# rate law this way removes every difference of rate-model convention between
# the two codes: what is left is the coupling itself.
#
# The oracle therefore needs no kinetics integrator at all — just one
# equilibrium solve per sample time, with the element totals set analytically.
#
#   conda run -n reaktoro-env python test/reference/reaktoro_coupling.py
#
# Paste the printed values into `REAKTORO_COUPLING` in test/coupling_reference.jl.

import reaktoro as rkt
import numpy as np

DB = "data/cemdata18-thermofun.json"

# Aqueous partition only: calcite is the kinetic species and is excluded.
SPECIES = "H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 CaOH+ Ca(CO3)@ Ca(HCO3)+".split()

N_H2O = 55.5          # mol
K_RATE = 1.0e-6       # mol/s, constant dissolution rate of calcite
TIMES = [0.0, 600.0, 1800.0, 3600.0]      # s

db = rkt.ThermoFunDatabase.fromFile(DB)
solution = rkt.AqueousPhase(" ".join(SPECIES))
solution.set(rkt.ActivityModelIdealAqueous())   # matches DiluteSolutionModel()
system = rkt.ChemicalSystem(db, solution)
names = [s.name() for s in system.species()]


def phi(xi):
    """Equilibrium speciation of the aqueous partition after `xi` mol of CaCO3
    has been released into it. Reaktoro conserves the element amounts of the
    initial state, so injecting Ca and C as Ca+2 + CO3-2 sets b_e exactly."""
    state = rkt.ChemicalState(system)
    state.temperature(25, "celsius")
    state.pressure(1, "bar")
    state.set("H2O@", N_H2O, "mol")
    if xi > 0.0:
        state.set("Ca+2", xi, "mol")
        state.set("CO3-2", xi, "mol")
    result = rkt.equilibrate(state)
    assert result.succeeded(), f"Reaktoro failed to equilibrate at xi={xi}"
    return np.array([state.speciesAmount(n) for n in names])


print(f"# Reaktoro {rkt.__version__}, Cemdata18, ideal activities")
print(f"# calcite dissolving at k = {K_RATE} mol/s")
print("const REAKTORO_COUPLING = Dict(")
for t in TIMES:
    n = phi(K_RATE * t)
    entries = ", ".join(f'"{nm}" => {n[k]:.9g}' for k, nm in enumerate(names))
    print(f"    {t} => Dict({entries}),")
print(")")

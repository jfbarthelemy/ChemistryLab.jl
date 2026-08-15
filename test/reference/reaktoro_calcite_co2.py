# Oracle for test/equilibrium_reference.jl.
#
# Reaktoro 2.13 reading the *same* Cemdata18 ThermoFun file this package ships,
# over the *same* species, under the *same* (ideal) activity model — the three
# knobs that have to match for a cross-code comparison to mean anything.
#
#   conda create -n reaktoro-env -c conda-forge reaktoro thermofun
#   conda run -n reaktoro-env python test/reference/reaktoro_calcite_co2.py
#
# Paste the printed values into `REAKTORO` in test/equilibrium_reference.jl.

import reaktoro as rkt
import numpy as np

DB = "data/cemdata18-thermofun.json"

SPECIES = "H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 CaOH+ Ca(CO3)@ Ca(HCO3)+".split()
N_H2O, N_CAL, N_CO2 = 55.5, 0.05, 0.01

db = rkt.ThermoFunDatabase.fromFile(DB)
solution = rkt.AqueousPhase(" ".join(SPECIES))
solution.set(rkt.ActivityModelIdealAqueous())   # matches DiluteSolutionModel()
system = rkt.ChemicalSystem(db, solution, rkt.MineralPhase("Cal"))
names = [s.name() for s in system.species()]


def speciate(n_co2):
    state = rkt.ChemicalState(system)
    state.temperature(25, "celsius")
    state.pressure(1, "bar")
    state.set("H2O@", N_H2O, "mol")
    state.set("Cal", N_CAL, "mol")
    state.set("CO2@", n_co2, "mol")
    result = rkt.equilibrate(state)
    assert result.succeeded(), "Reaktoro failed to equilibrate"
    return np.array([state.speciesAmount(n) for n in names])


amounts = speciate(N_CO2)

# Central differences, and the spread across step sizes that bounds how well
# this oracle can be trusted in the first place.
derivs = {}
for h in (1e-3, 1e-4, 1e-5):
    derivs[h] = (speciate(N_CO2 + h) - speciate(N_CO2 - h)) / (2 * h)
spread = np.max(np.abs(derivs[1e-3] - derivs[1e-5]))
d = derivs[1e-5]

print(f"# Reaktoro {rkt.__version__}, Cemdata18, ideal activities")
print(f"# finite-difference spread across h: {spread:.3g}")
print("const REAKTORO = Dict(")
for k, name in enumerate(names):
    print(f'    "{name}" => (n = {amounts[k]:.9g}, dn = {d[k]:.9g}),')
print(")")

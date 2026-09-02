# [Effect of Water/Cement Ratio on Cement Hydration](@id sec-wc-ratio)

The **water-to-cement ratio** (w/c) is the single most important mix-design parameter
of concrete. It controls workability, compressive strength, and durability simultaneously.
From a thermodynamic perspective, it determines how much water is available to hydrate the
clinker phases, which in turn governs the nature and amount of hydration products, the
porosity of the paste, and the pH of the pore solution.

This example scans w/c from 0.30 to 0.60 and tracks, at **full thermodynamic
equilibrium**, the pH of the pore solution, the hydrate assemblage, and the
porosity of the hardened paste in the sealed-curing convention.

!!! warning "Equilibrium answers a narrower question than mix design does"
    Read the scan below for what it is: the assemblage a paste would reach if
    every reaction ran to completion. It is **not** the state of a real paste at
    an age, and three of the effects usually attributed to w/c are absent from it
    by construction — no clinker survives at any w/c, there is no water-limited
    regime, and there is no optimum. Those are **kinetic** limitations, carried
    by [`powers_alpha_max`](@ref) and the rate laws of the
    [kinetics tutorial](@ref sec-kinetics), not by the Gibbs minimum. This page
    is the reference the kinetic calculation converges toward; the
    [coupled hydration example](@ref sec-coupled-hydration) is the one to read
    for an age.

---

## System setup

The same clinker composition and species set as the "simplified clinker dissolution" example are used here.

```@example wc_setup
using ChemistryLab
using DynamicQuantities

substances = build_species(datapath("cemdata18-thermofun.json"))

input_species = split("C3S C2S C3A C4AF Gp Anh Portlandite Jennite H2O@ ettringite monosulphate12 C3AH6 C3FH6 C4FH13")
species = speciation(substances, input_species; aggregate_state = [AS_AQUEOUS])

cs = ChemicalSystem(species, CEMDATA_PRIMARIES)
```

The clinker composition (mass fractions of the anhydrous cement phases) is fixed throughout the scan:

| Phase | Symbol | Mass fraction |
|:------|:-------|:--------------|
| Alite | `C3S`  | 67.8 % |
| Belite | `C2S` | 16.6 % |
| Aluminate | `C3A` | 4.0 % |
| Ferrite | `C4AF` | 7.2 % |
| Gypsum | `Gp`  | 2.8 % |

```@example wc_setup
compo = ["C3S" => 0.678, "C2S" => 0.166, "C3A" => 0.040, "C4AF" => 0.072, "Gp" => 0.028]
c     = sum(last.(compo))   # cement mass fraction (= 0.984 here)
```

---

## Building the equilibrium solver

A single [`EquilibriumSolver`](@ref) is compiled once and reused for every w/c point:

```@example wc_setup
using Optimization, OptimizationIpopt

opt = IpoptOptimizer(
    acceptable_tol        = 1e-10,
    dual_inf_tol          = 1e-10,
    acceptable_iter       = 100,
    constr_viol_tol       = 1e-10,
    warm_start_init_point = "no",
)

solver = EquilibriumSolver(
    cs,
    DiluteSolutionModel(),
    opt;
    variable_space = Val(:linear),
    abstol  = 1e-8,
    reltol  = 1e-8,
)
```

---

## Scanning the w/c ratio

For each value of w/c the fresh state is rebuilt from scratch, and the total mass is
normalized to 1 kg of paste (cement + water) so that all amounts are comparable
across the scan. The fresh state is kept: it is the volume reference the porosity
is referred to.

```@example wc_setup
sp_idx   = Dict(symbol(s) => i for (i, s) in enumerate(cs.species))
wc_range = range(0.30, 0.60; length = 13)

function fresh_paste(wc)
    w, mtot = wc * c, c + wc * c
    st = ChemicalState(cs)
    for (sym, mfrac) in compo
        set_quantity!(st, sym, mfrac / mtot * u"kg")
    end
    set_quantity!(st, "H2O@", w / mtot * u"kg")
    V = volume(st)
    set_quantity!(st, "H+",  1e-7u"mol/L" * V.liquid)   # charge seed, pH-neutral
    set_quantity!(st, "OH-", 1e-7u"mol/L" * V.liquid)
    return st
end

pH_vals    = Float64[]
ϕ_liquid   = Float64[]
ϕ_void     = Float64[]
ϕ_total    = Float64[]
n_portl    = Float64[]
n_mono     = Float64[]
n_ett      = Float64[]
n_jennite  = Float64[]
n_clinker  = Float64[]

for wc in wc_range
    fresh = fresh_paste(wc)
    eq    = solve(solver, deepcopy(fresh))   # solve mutates its argument
    ϕ     = porosity(eq, fresh)              # sealed-curing convention
    amount(sym) = ustrip(eq.n[sp_idx[sym]])
    push!(pH_vals,   pH(eq))
    push!(ϕ_liquid,  ϕ.liquid)
    push!(ϕ_void,    ϕ.void)
    push!(ϕ_total,   ϕ.total)
    push!(n_portl,   amount("Portlandite"))
    push!(n_mono,    amount("monosulphate12"))
    push!(n_ett,     amount("ettringite"))
    push!(n_jennite, amount("Jennite"))
    push!(n_clinker, sum(amount(s) for s in ("C3S", "C2S", "C3A", "C4AF")))
end
nothing # hide
```

!!! note "Why `porosity(eq, fresh)` and not `porosity(eq)`"
    The one-argument [`porosity`](@ref) divides the pore volume by the **current**
    total volume, which is right for a fixed-volume aqueous system and wrong for a
    setting binder: hydration products occupy less space than the reactants they
    consume, so the paste's own volume shrinks under the ratio. The two-argument
    form refers everything to the fresh volume and splits the result into the
    water-filled porosity and the empty porosity that the Le Chatelier
    contraction creates. On this mix at w/c = 0.50 the one-argument form returns
    0.285 against a total of 0.338, a gap of 5.3 points of porosity, the total
    volume having shrunk by 7.4 %.

---

## Results

### pH and porosity

```@example wc_setup
using Plots

p1 = plot(collect(wc_range), pH_vals;
    xlabel = "w/c ratio", ylabel = "Pore solution pH", label = "pH",
    linewidth = 2, marker = :circle, markersize = 4, color = :steelblue,
    title = "Pore solution pH", ylims = (11.5, 13.5), legend = :bottomright)

p2 = plot(collect(wc_range), ϕ_total .* 100;
    xlabel = "w/c ratio", ylabel = "Porosity (%)", label = "total",
    linewidth = 2, marker = :circle, markersize = 4, color = :firebrick,
    title = "Porosity referred to the fresh volume",
    ylims = (0, 50), legend = :topleft)
plot!(p2, collect(wc_range), ϕ_liquid .* 100;
    label = "water-filled", linewidth = 2, marker = :square, markersize = 3,
    color = :steelblue)
plot!(p2, collect(wc_range), ϕ_void .* 100;
    label = "empty (Le Chatelier)", linewidth = 2, marker = :diamond,
    markersize = 3, color = :seagreen)

plot(p1, p2; layout = (1, 2), left_margin = 8Plots.mm,
     bottom_margin = 8Plots.mm, size = (950, 410))
```

The pH is **flat to three decimals** over the whole range. That is the signature
of a buffered solution: portlandite is present at every w/c, so the calcium and
hydroxide activities are pinned by its saturation, and in a dilute-solution model
those activities do not know how large the pore volume is. Diluting a saturated
solution with more of its own solvent does not change its pH — it dissolves more
portlandite. The pH would start to move only once portlandite is exhausted, which
this composition never does.

### Phase assemblage

```@example wc_setup
p3 = plot(collect(wc_range), n_portl;
    xlabel = "w/c ratio", ylabel = "Amount (mol / kg of paste)",
    label = "Portlandite  Ca(OH)₂", linewidth = 2, marker = :circle,
    markersize = 4, color = :steelblue, title = "Hydrate assemblage at equilibrium",
    legend = :right)
plot!(p3, collect(wc_range), n_jennite;
    label = "Jennite (C-S-H)", linewidth = 2, marker = :square,
    markersize = 3, color = :firebrick)
plot!(p3, collect(wc_range), n_mono;
    label = "monosulphate12 (AFm)", linewidth = 2, marker = :diamond,
    markersize = 3, color = :seagreen)
plot!(p3, collect(wc_range), n_ett;
    label = "ettringite (AFt)", linewidth = 2, marker = :utriangle,
    markersize = 3, color = :darkorange)
plot(p3; left_margin = 8Plots.mm, bottom_margin = 8Plots.mm, size = (700, 420))
```

```@example wc_setup
using Printf
@printf "ettringite, largest value over the scan : %.3e mol\n" maximum(n_ett)
@printf "clinker left, largest value over the scan : %.3e mol\n" maximum(n_clinker)
```

---

## Analysis

| Quantity | What the scan gives | Why |
|:--|:--|:--|
| **pH** | 12.39, constant to three decimals | buffered by portlandite saturation; independent of pore volume in a dilute model |
| **Portlandite, C-S-H, AFm** | decrease with w/c, in proportion to the cement fraction | amounts are per kg of *paste*; adding water dilutes the binder, it does not change what a gram of cement produces |
| **Ettringite** | at the solver's lower bound, ``10^{-8}`` mol, i.e. absent | see below |
| **Clinker left** | at the same bound at every w/c, 0.30 included | see below |
| **Total porosity** | 12.2 % to 41.0 %, monotone, no optimum | the excess water has nowhere to go but the pore space |
| **Empty porosity** | 9.8 % down to 6.6 % | the Le Chatelier contraction is roughly fixed per gram of cement, so it is a smaller fraction of a larger reference volume |

Two of those deserve stating plainly, because the intuition of mix design points
the other way.

!!! warning "Ettringite does not form here, and that is correct"
    AFt needs about three sulfates per aluminate; this clinker has 2.8 % gypsum
    against 4.0 % C₃A plus the ferrite, so at **equilibrium** the sulfate is all
    taken up by the AFm phase `monosulphate12`, and `ettringite` stays at the
    solver's lower bound — ``1.3 \times 10^{-8}`` mol at most over the whole
    scan, which is absence, not a small amount. Ettringite is the phase that forms **early**, while
    sulfate is still locally abundant, and then converts to AFm as it runs out.
    It is a kinetic intermediate for this mix, so an equilibrium scan cannot show
    it. Raise the gypsum content and AFt becomes stable — that is the
    sulfate-balance calculation, and it is worth doing before reading anything
    into an AFt amount.

!!! warning "No clinker survives, so there is no water-limited regime"
    A Gibbs minimization consumes the clinker at every w/c, including 0.30 — the
    four anhydrous phases together never exceed ``1.7 \times 10^{-8}`` mol, again
    the lower bound. The
    reason is that the anhydrous phases are never the stable ones in the presence
    of water, and the minimum can always form a less hydrous assemblage rather
    than leave alite standing. The **water-limited regime** of mix design — the
    unreacted clinker of a low-w/c paste, and the maximum degree of hydration
    ``\alpha_{\max} \le w/c \,/\, 0.42`` of [Powers1948](@cite) — is a
    statement about **rates and access to water**, not about the Gibbs minimum.
    It enters ChemistryLab through [`powers_alpha_max`](@ref) in the kinetic rate
    laws. Consequently this scan shows **no optimum w/c and no inflection**: the
    porosity rises monotonically, and the minimum-porosity mix design does not
    appear. Looking for it here is looking in the wrong calculation.

### Assumptions behind these numbers

  - **Complete reaction.** Full equilibrium, no time, no kinetic barrier.
  - **Sealed curing.** No water exchanged with the outside, so the empty
    porosity from the chemical shrinkage stays empty. An immersed specimen would
    draw water in and `ϕ.void` would fill.
  - **The fresh volume is the reference**, and it is held fixed — the specimen
    keeps its cast dimensions, the contraction showing up as internal void rather
    than as shrinkage of the outside. This is the usual convention for a set
    paste; it is wrong before setting, when the material still contracts
    externally.
  - **Ideal molar volumes.** Phase volumes are the sum of ``n_i V_i^0``, with no
    mixing term.
  - **Dilute solution model**, so no ionic-strength correction. The pore solution
    of a cement paste is around 0.1–0.3 mol/kg, where activity coefficients
    depart from unity by tens of percent — use [`HKFActivityModel`](@ref) or
    [`DaviesActivityModel`](@ref) if the ion concentrations themselves matter.
    The pH being buffered, it is the quantity least affected by this choice.
  - **The species list is closed.** Only the 14 species selected above may form;
    siliceous hydrogarnet, hydrotalcite and the alkali sulfates are absent, as
    are the alkalis themselves, which in a real paste raise the pore-solution pH
    to 13 or above.

!!! note "Extending the scan"
    To study **supplementary cementitious materials**, substitute part of the
    clinker and add the corresponding species from `cemdata18` (e.g. `C2ASH8`).
    For **carbonation**, add `CO2@`, `HCO3-`, `CO3-2` and the carbonate phases —
    see the [cement carbonation example](@ref sec-cement-carbonation).

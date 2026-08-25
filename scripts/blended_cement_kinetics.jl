# =============================================================================
# blended_cement_kinetics.jl
#
# Kinetic simulation of blended cement hydration with mineral additions:
#   - 63% OPC clinker (CEM I 52.5 R, 4 Parrot-Killoh phases)
#   - 30% ground granulated blast-furnace slag (GGBS, slow)
#   - 7% metakaolin (MK, more reactive than slag)
#
# CEMDATA18 does not contain GGBS or MK as reactants: we create custom Species
# with a dummy ΔₐG⁰ (the Parrot-Killoh model ignores Ω).
# The reference mole of each addition is defined by the explicit :M field
# of the Species (mass per "representative formula unit"), which must be
# consistent with ΔᵣH⁰ for correct calorimetry:
#
#     |ΔᵣH⁰| [J/mol] = specific_heat [J/g] × M [g/mol]
#
# Usage:
#   julia --project=scripts scripts/blended_cement_kinetics.jl
#   or from the REPL:  include("scripts/blended_cement_kinetics.jl")
#   or, in VS Code, just run the file: it activates `scripts/` itself, and every
#   data file it reads is named with `datapath`, so no working directory matters.

import Pkg
Pkg.activate(@__DIR__; io = devnull)

using ChemistryLab
using OrdinaryDiffEq
using DynamicQuantities
using OrderedCollections
using Printf

# ── 1. Base ChemicalSystem from CEMDATA18 ────────────────────────────────────
#
# We select the clinker phases, the main hydration products
# (portlandite, C-S-H, ettringite, AFm, hydrotalcite, stratlingite...) and water.
# GGBS and MK will be added afterwards as custom species.

const DATA_FILE = datapath("cemdata18-thermofun.json")

substances = build_species(DATA_FILE)

input_species = split(
    "C3S C2S C3A C4AF " *
        "Portlandite Jennite ettringite monosulphate12 C3AH6 C3FH6 " *
        "straetlingite hydrotalcite " *
        "H2O@",
)

species_base = speciation(substances, input_species; aggregate_state = [AS_AQUEOUS])
cs_base = ChemicalSystem(species_base, CEMDATA_PRIMARIES)

@info "Base system: $(length(cs_base.species)) species"

# ── 2. Custom Species for GGBS and MK ────────────────────────────────────────
#
# GGBS (ground granulated blast-furnace slag):
#   Typical CEM I-like composition: ~42% CaO, 35% SiO₂, 12% Al₂O₃, 8% MgO.
#   Representative formula (charge-neutral): CaAl₂Si₂O₈ (anorthite analogy).
#   Molar mass of the formula unit overridden to 95 g/mol so that:
#       |ΔᵣH⁰| = 380 J/g × 95 g/mol ≈ 36 100 J/mol  (Gruyaert 2010)
#
# MK (metakaolin):
#   Exact formula: Al₂Si₂O₇ (dehydroxylated kaolinite).
#   Auto-computed molar mass: M = 222 g/mol.
#   |ΔᵣH⁰| = 250 J/g × 222 g/mol ≈ 55 500 J/mol  (Lothenbach 2011)
#
# Dummy ΔₐG⁰: very negative → Ω ≈ 0 (dissolution always favored).
# Parrot-Killoh ignores Ω; the value does not affect kinetic rates.

const _dummy_G = NumericFunc((T, P) -> -1_200_000.0, (:T, :P), u"J/mol")

# GGBS: formula CaAl₂Si₂O₈, :M overridden to 95 g/mol
sp_ggbs = Species(
    "CaAl2Si2O8";
    symbol = "GGBS",
    name = "GGBS",
    aggregate_state = AS_CRYSTAL,
    properties = Dict{Symbol, Any}(
        :M => 0.095u"kg/mol",
        :ΔₐG⁰ => _dummy_G,
    ),
)

# MK: formula Al₂Si₂O₇ (metakaolin), M = 222 g/mol auto-computed
sp_mk = Species(
    "Al2Si2O7";
    symbol = "MK",
    name = "MK",
    aggregate_state = AS_CRYSTAL,
    properties = Dict{Symbol, Any}(
        :ΔₐG⁰ => _dummy_G,
    ),
)

# ── 3. Extended ChemicalSystem ───────────────────────────────────────────────

all_species = vcat(cs_base.species, sp_ggbs, sp_mk)
cs = ChemicalSystem(all_species, CEMDATA_PRIMARIES)

@info "Extended system: $(length(cs.species)) species (including GGBS and MK)"

# ── 4. Initial state ─────────────────────────────────────────────────────────
#
# Ternary cement (1 kg): 63% OPC clinker, 30% GGBS slag, 7% metakaolin
# Water/binder ratio w/b = 0.40
#
# Clinker phases (mass fractions in CEM I 52.5 R clinker):
#   C₃S: 61.9%, C₂S: 16.5%, C₃A: 8.0%, C₄AF: 8.7%   (Lavergne 2018)
# Clinker fraction in the ternary cement: 0.63

const WB = 0.4      # water/binder ratio

const CLINKER_FRAC = 0.63
const GGBS_FRAC = 0.3
const MK_FRAC = 0.07

const CLINKER_COMP = (C3S = 0.619, C2S = 0.165, C3A = 0.08, C4AF = 0.087)

const COMPOSITION = (
    C3S = CLINKER_COMP.C3S * CLINKER_FRAC,
    C2S = CLINKER_COMP.C2S * CLINKER_FRAC,
    C3A = CLINKER_COMP.C3A * CLINKER_FRAC,
    C4AF = CLINKER_COMP.C4AF * CLINKER_FRAC,
    GGBS = GGBS_FRAC,
    MK = MK_FRAC,
)

state0 = ChemicalState(cs)
for (name, frac) in pairs(COMPOSITION)
    set_quantity!(state0, string(name), frac * u"kg")
end
set_quantity!(state0, "H2O@", WB * u"kg")

# ── 5. Parrot & Killoh kinetic models ────────────────────────────────────────
#
# α_max according to Powers (1948): hydration limited by available water
const α_max = min(1.0, WB / 0.42)

# Clinker (Parrot & Killoh 1984, Schindler & Folliard 2005 corrections)
pk_C3S = parrot_killoh(PK_PARAMS_C3S, "C3S"; α_max)
pk_C2S = parrot_killoh(PK_PARAMS_C2S, "C2S"; α_max)
pk_C3A = parrot_killoh(PK_PARAMS_C3A, "C3A"; α_max)
pk_C4AF = parrot_killoh(PK_PARAMS_C4AF, "C4AF"; α_max)

# GGBS (slag): PK parameters adapted from the literature
#   Low K₁ → slow initial kinetics (vitreous, little nucleation)
#   High Ea → strong temperature sensitivity
#   α_max = 0.90 (vitreous, partial hydration)
#   References: Richardson & Groves (1992), Chen & Brouwers (2007)
pk_ggbs = parrot_killoh(
    (
        K₁ = 0.15u"1/d", N₁ = 2.0, K₂ = 0.003u"1/d", N₂ = 2.0,
        K₃ = 0.0015u"1/d", N₃ = 3.5, B = 0.2,
        Ea = 46_000.0u"J/mol", T_ref = 293.15u"K",
    ),
    "GGBS";
    α_max = 0.9,
)

# MK (metakaolin): PK parameters adapted from the literature
#   High K₁ → fast initial reactivity (high specific surface area)
#   α_max = 0.95 (nearly complete under normal conditions)
#   References: Lothenbach et al. (2011), Deschner et al. (2012)
pk_mk = parrot_killoh(
    (
        K₁ = 0.7u"1/d", N₁ = 1.5, K₂ = 0.008u"1/d", N₂ = 2.0,
        K₃ = 0.002u"1/d", N₃ = 3.5, B = 0.3,
        Ea = 48_000.0u"J/mol", T_ref = 293.15u"K",
    ),
    "MK";
    α_max = 0.95,
)

# ── 6. Kinetic reaction list ─────────────────────────────────────────────────
#
# Balanced clinker reactions (CEMDATA18).
# Jennite in CEMDATA18 = (SiO₂)(CaO)_{5/3}(H₂O)_{21/10}, Ca:Si = 5/3.
#
#   C₃S  + 103/30 H₂O  →  Jennite  + 4/3 Portlandite       (ΔᵣH⁰ ≈ −124 kJ/mol)
#   C₂S  +  73/30 H₂O  →  Jennite  + 1/3 Portlandite       (ΔᵣH⁰ ≈  −48 kJ/mol)
#   C₃A  +  6     H₂O  →  C₃AH₆                            (ΔᵣH⁰ ≈ −261 kJ/mol)
#   C₄AF + 2 Portlandite + 10 H₂O → C₃AH₆ + C₃FH₆         (ΔᵣH⁰ ≈ −147 kJ/mol)
#
# GGBS and MK have artificial formulas (no ΔₐH⁰ in database); ΔᵣH⁰ is set
# directly on the reaction (thermodynamic convention: negative = exothermic).
#   GGBS: 380 J/g ×  95 g/mol ≈  36 100 J/mol  (Gruyaert 2010)
#   MK  : 250 J/g × 222 g/mol ≈  55 500 J/mol  (Lothenbach 2011)

sp(name) = cs[name]

rxn_C3S = Reaction(
    OrderedDict(sp("C3S") => 1.0, sp("H2O@") => 103 / 30),
    OrderedDict(sp("Jennite") => 1.0, sp("Portlandite") => 4 / 3);
    symbol = "C₃S hydration",
)
rxn_C3S[:rate] = pk_C3S

rxn_C2S = Reaction(
    OrderedDict(sp("C2S") => 1.0, sp("H2O@") => 73 / 30),
    OrderedDict(sp("Jennite") => 1.0, sp("Portlandite") => 1 / 3);
    symbol = "C₂S hydration",
)
rxn_C2S[:rate] = pk_C2S

rxn_C3A = Reaction(
    OrderedDict(sp("C3A") => 1.0, sp("H2O@") => 6.0),
    OrderedDict(sp("C3AH6") => 1.0);
    symbol = "C₃A hydration",
)
rxn_C3A[:rate] = pk_C3A

rxn_C4AF = Reaction(
    OrderedDict(sp("C4AF") => 1.0, sp("Portlandite") => 2.0, sp("H2O@") => 10.0),
    OrderedDict(sp("C3AH6") => 1.0, sp("C3FH6") => 1.0);
    symbol = "C₄AF hydration",
)
rxn_C4AF[:rate] = pk_C4AF

# GGBS: artificial formula → approximate products (Jennite + stratlingite)
# ΔᵣH⁰ set directly (no ΔₐH⁰ for custom species)
rxn_GGBS = Reaction(
    OrderedDict(sp("GGBS") => 1.0, sp("H2O@") => 3.0),
    OrderedDict(sp("Jennite") => 0.05, sp("straetlingite") => 0.2);
    symbol = "GGBS hydration",
)
rxn_GGBS[:rate] = pk_ggbs
rxn_GGBS[:ΔᵣH⁰] = NumericFunc((T) -> -36_100.0, (:T,), u"J/mol")

# MK: Al₂Si₂O₇ + 2 Ca(OH)₂ + 5 H₂O → stratlingite
rxn_MK = Reaction(
    OrderedDict(sp("MK") => 1.0, sp("Portlandite") => 2.0, sp("H2O@") => 5.0),
    OrderedDict(sp("straetlingite") => 1.0);
    symbol = "MK hydration",
)
rxn_MK[:rate] = pk_mk
rxn_MK[:ΔᵣH⁰] = NumericFunc((T) -> -55_500.0, (:T,), u"J/mol")

kinetic_reactions = [rxn_C3S, rxn_C2S, rxn_C3A, rxn_C4AF, rxn_GGBS, rxn_MK]

# ── 7. Kinetics problem ─────────────────────────────────────────────────────
#
# tspan: 28 days → shows the late reactivity of slag
# equilibrium_solver = nothing: no aqueous re-speciation (sufficient for PK)

const TSPAN = (0.0u"s", 28.0u"d")

# ── 8. Semi-adiabatic calorimeter ────────────────────────────────────────────
#
# Cp [J/K]: ternary cement + water + Dewar flask
#   cement: 1 kg × 800 J/(kg·K)
#   water : WB kg × 4186 J/(kg·K)
#   Dewar : 1 kg × 900 J/(kg·K)
# Quadratic losses (Lavergne et al. 2018)

cal = SemiAdiabaticCalorimeter(;
    Cp = (1.0 * 800.0 + WB * 4186.0 + 1.0 * 900.0) * u"J/K",
    T_env = 293.15u"K",
    heat_loss = ΔT -> 0.3 * ΔT + 0.003 * ΔT^2,
    T0 = 293.15u"K",
)

kp = KineticsProblem(
    cs,
    kinetic_reactions,
    state0,
    TSPAN;
    calorimeter = cal,
    equilibrium_solver = nothing,
)

# ── 9. Integration ───────────────────────────────────────────────────────────

@info "Integration in progress (28 days, Rodas5P)..."
ks = KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-6, abstol = 1.0e-9)
sol = integrate(kp, ks)
@info "Done: $(length(sol.t)) accepted steps."

# ── 10. Post-processing ─────────────────────────────────────────────────────

t_T, T_K_vec = temperature_profile(sol, cal)
t_Q, Q_J_vec = cumulative_heat(sol, cal)
T_°C_vec = T_K_vec .- 273.15
Q_kJ_vec = Q_J_vec ./ 1000.0

n0_kin = [sol.prob.p.n_initial_full[i] for i in kp.idx_kinetic]
n_kin = [[u[i] for u in sol.u] for i in eachindex(n0_kin)]
t_h = sol.t ./ 3600.0

function phase_alpha(cs, kp, sol, n0_kin, n_kin, name)
    sp_idx = findfirst(sp -> ChemistryLab.symbol(sp) == name, cs.species)
    pos = findfirst(==(sp_idx), kp.idx_kinetic)
    isnothing(pos) && return fill(NaN, length(sol.t))
    return 1.0 .- n_kin[pos] ./ n0_kin[pos]
end

α_C3S = phase_alpha(cs, kp, sol, n0_kin, n_kin, "C3S")
α_C2S = phase_alpha(cs, kp, sol, n0_kin, n_kin, "C2S")
α_C3A = phase_alpha(cs, kp, sol, n0_kin, n_kin, "C3A")
α_C4AF = phase_alpha(cs, kp, sol, n0_kin, n_kin, "C4AF")
α_GGBS = phase_alpha(cs, kp, sol, n0_kin, n_kin, "GGBS")
α_MK = phase_alpha(cs, kp, sol, n0_kin, n_kin, "MK")

# Mass-weighted mean degree of hydration
w = COMPOSITION
α_mean = (
    w.C3S .* α_C3S .+ w.C2S .* α_C2S .+
        w.C3A .* α_C3A .+ w.C4AF .* α_C4AF .+
        w.GGBS .* α_GGBS .+ w.MK .* α_MK
) ./ (w.C3S + w.C2S + w.C3A + w.C4AF + w.GGBS + w.MK)

# ── Summary ──────────────────────────────────────────────────────────────────

idx_7d = findlast(t -> t <= 7 * 86400, sol.t)

println()
println("╔═══════════════════════════════════════════════════════════╗")
println("║   Ternary cement kinetics — ChemistryLab                  ║")
@printf "║   63%% CEM I + 30%% GGBS + 7%% MK  │  w/b = %.2f  │  28 d   ║\n" WB
println("╠═══════════════════════════════════════════════════════════╣")
@printf "║  ΔT max      = %6.2f °C                                  ║\n" maximum(T_°C_vec) - 20.0
println("╠═══════════════════════════════════════════════════════════╣")
@printf "║  α(C₃S)  7d  = %6.4f   28d = %6.4f                    ║\n" α_C3S[idx_7d] α_C3S[end]
@printf "║  α(C₂S)  7d  = %6.4f   28d = %6.4f                    ║\n" α_C2S[idx_7d] α_C2S[end]
@printf "║  α(GGBS) 7d  = %6.4f   28d = %6.4f                    ║\n" α_GGBS[idx_7d] α_GGBS[end]
@printf "║  α(MK)   7d  = %6.4f   28d = %6.4f                    ║\n" α_MK[idx_7d] α_MK[end]
@printf "║  ᾱ mean   7d  = %6.4f   28d = %6.4f                    ║\n" α_mean[idx_7d] α_mean[end]
println("╠═══════════════════════════════════════════════════════════╣")
@printf "║  Q  7d = %6.1f kJ/kg   Q 28d = %6.1f kJ/kg            ║\n" Q_kJ_vec[idx_7d] Q_kJ_vec[end]
println("╚═══════════════════════════════════════════════════════════╝")

# ── 11. Plots ────────────────────────────────────────────────────────────────

using Plots
gr()

p1 = plot(
    t_T ./ 3600, T_°C_vec;
    xlabel = "Time [h]", ylabel = "T [°C]",
    title = "Temperature (semi-adiabatic calorimeter)",
    label = "T(t)", lw = 2, color = :red,
)
hline!(p1, [20.0]; linestyle = :dash, color = :gray, label = "T₀ = T_env")

p2 = plot(
    t_h, [α_C3S α_C2S α_C3A α_C4AF α_GGBS α_MK α_mean];
    xlabel = "Time [h]", ylabel = "Degree of hydration α",
    title = "Phase hydration",
    label = ["C₃S" "C₂S" "C₃A" "C₄AF" "GGBS" "MK" "ᾱ mean"],
    lw = 2, ls = [:solid :dash :dot :dashdot :solid :dash :solid],
    color = [:blue :cyan :green :orange :brown :purple :black],
)
hline!(p2, [α_max]; linestyle = :dash, color = :black, label = "α_max (Powers)")

p3 = plot(
    t_Q ./ 3600, Q_kJ_vec;
    xlabel = "Time [h]", ylabel = "Q [kJ/kg cement]",
    title = "Cumulative heat",
    label = "Q(t)", lw = 2, color = :purple,
)

display(
    plot(
        p1, p2, p3;
        layout = (1, 3), top_margin = 7Plots.mm, left_margin = 8Plots.mm, bottom_margin = 8Plots.mm, size = (1500, 450),
        plot_title = "ChemistryLab — Ternary cement 63% OPC + 30% GGBS + 7% MK  (w/b=$WB)",
    )
)

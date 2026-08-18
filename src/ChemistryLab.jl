# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

__precompile__(true)

"""
    ChemistryLab

Top-level module for parsing, representing and manipulating chemical
formulas, species, stoichiometric matrices and ThermoFun / PHREEQC-like data.

# Overview

  - Canonical containers: `AtomGroup`, `Formula`.
  - Species representations: `Species`, `CemSpecies`.
  - Parsers and IO helpers for ThermoFun / PHREEQC data (substances, reactions).
  - Stoichiometric matrix construction and reaction helpers.
  - Utilities for units, colored and Unicode terminal output.

# Examples

```julia
julia> using ChemistryLab

julia> f = Formula("H2O")
Formula{Int64}
    formula: H2O ◆ H₂O
composition: H => 2, O => 1
     charge: 0

julia> f[:H]
2

julia> H2O = Species(f; name="Vapour", symbol="H₂O", aggregate_state=AS_GAS, class=SC_GASFLUID)
Species{Int64}
           name: Vapour
         symbol: H₂O
        formula: H2O ◆ H₂O
          atoms: H => 2, O => 1
         charge: 0
aggregate_state: AS_GAS
          class: SC_GASFLUID
     properties: M = 0.0180149999937744 kg mol⁻¹

julia> CO2 = Species(
           Dict(:C=>1, :O=>2);
           name="Carbon dioxide",
           symbol="CO₂",
           aggregate_state=AS_GAS,
           class=SC_GASFLUID,
       ) # definition from Dict
Species{Int64}
           name: Carbon dioxide
         symbol: CO₂
        formula: CO2 ◆ CO₂
          atoms: O => 2, C => 1
         charge: 0
aggregate_state: AS_GAS
          class: SC_GASFLUID
     properties: M = 0.04400899998479143 kg mol⁻¹

julia> O2 = Species("O2"; name="Dioxygen", symbol="O₂", aggregate_state=AS_GAS, class=SC_GASFLUID) # definition from String
Species{Int64}
           name: Dioxygen
         symbol: O₂
        formula: O2 ◆ O₂
          atoms: O => 2
         charge: 0
aggregate_state: AS_GAS
          class: SC_GASFLUID
     properties: M = 0.03199799998894218 kg mol⁻¹

julia> C3H8 = Species(
           "C₃H₈"; name="Propane", symbol="C₃H₈", aggregate_state=AS_GAS, class=SC_GASFLUID
       ) # definition from Unicode String
Species{Int64}
           name: Propane
         symbol: C₃H₈
        formula: C₃H₈ ◆ C3H8
          atoms: C => 3, H => 8
         charge: 0
aggregate_state: AS_GAS
          class: SC_GASFLUID
     properties: M = 0.04409699998476102 kg mol⁻¹

julia> r = Reaction([C3H8, O2, CO2, H2O])
  equation: C₃H₈ + 5O₂ = 4H₂O + 3CO₂
 reactants: C₃H₈ => 1, O₂ => 5
  products: H₂O => 4, CO₂ => 3
    charge: 0
```
"""
module ChemistryLab

    import Base: ==, +, -, *, /, //, ^

    using RuntimeGeneratedFunctions

    RuntimeGeneratedFunctions.init(@__MODULE__)

    include("utils/misc.jl")
    include("utils/subsuperscripts.jl")

    include("thermodynamics/thermo_factories.jl")
    include("thermodynamics/water_properties.jl")
    include("thermodynamics/thermo_models.jl")
    include("thermodynamics/precompile.jl")

    include("chemical_structs/element_order.jl")
    include("chemical_structs/parsing_tools.jl")
    include("chemical_structs/formulas.jl")
    include("chemical_structs/species.jl")
    include("chemical_structs/solid_solutions.jl")
    include("chemical_structs/reactions.jl")
    include("chemical_structs/speciation.jl")
    include("chemical_structs/stoich_matrices.jl")
    include("chemical_structs/chemical_systems.jl")
    include("chemical_structs/chemical_states.jl")
    include("chemical_structs/volume_fractions.jl")

    include("databases/phreeqc_dat.jl")
    include("databases/thermofun_json.jl")
    include("databases/merge_dat_json.jl")

    include("equilibrium/activities.jl")
    include("equilibrium/equilibrium_problems.jl")
    include("equilibrium/equilibrium_solver.jl")
    include("equilibrium/dual_solver.jl")

    include("kinetics/rate_models.jl")
    include("kinetics/kinetics_reactions.jl")
    include("kinetics/kinetics_problems.jl")
    include("kinetics/kinetics_solver.jl")
    include("kinetics/calorimetry.jl")
    include("kinetics/kinetics_postprocessing.jl")

    export SymbolicFunc,
        ThermoFactory,
        NumericFunc,
        AbstractFunc,
        derivative

    export WaterThermoProps,
        WaterElectroProps,
        HKFGState,
        SpeciesElectroPropsHKF,
        water_thermo_props,
        water_electro_props_jn,
        hkf_g_function,
        species_electro_props_hkf

    export THERMO_MODELS,
        THERMO_FACTORIES,
        add_thermo_model,
        build_thermo_functions

    export ATOMIC_ORDER,
        CEMENT_TO_MENDELEEV,
        OXIDE_ORDER,
        CEMDATA_PRIMARIES

    export stoich_coef_round,
        phreeqc_to_unicode,
        unicode_to_phreeqc,
        colored_formula,
        parse_formula,
        extract_charge,
        to_mendeleev,
        parse_equation,
        colored_equation,
        format_equation

    export AtomGroup,
        Formula,
        expr,
        phreeqc,
        unicode,
        colored,
        composition,
        charge,
        check_mendeleev,
        calculate_molar_mass,
        stoichtype,
        pprint

    export AggregateState,
        AS_UNDEF,
        AS_AQUEOUS,
        AS_CRYSTAL,
        AS_GAS

    export Class,
        SC_UNDEF,
        SC_AQSOLVENT,
        SC_AQSOLUTE,
        SC_COMPONENT,
        SC_GASFLUID,
        SC_SSENDMEMBER

    export AbstractSpecies,
        Species,
        CemSpecies,
        name,
        symbol,
        formula,
        cemformula,
        atoms,
        atoms_charge,
        oxides,
        oxides_charge,
        components,
        aggregate_state,
        class,
        properties,
        apply

    export AbstractReaction,
        Reaction,
        CemReaction,
        reactants,
        products,
        charge,
        equation,
        equal_sign,
        simplify_reaction
    @eval export $(Symbol.(EQUAL_OPS)...)

    export union_atoms, speciation

    export StoichMatrix,
        CanonicalStoichMatrix,
        pull_primaries,
        push_primaries,
        mass_matrix,
        reactions

    export DualEquilibriumSolver,
        optimality_certificate,
        AbstractSolidSolutionModel,
        IdealSolidSolutionModel,
        RedlichKisterModel,
        AbstractSolidSolutionPhase,
        SolidSolutionPhase,
        end_members,
        model,
        with_class

    export ChemicalSystem,
        aqueous,
        crystal,
        gas,
        solutes,
        solvent,
        components,
        gasfluid,
        get_reaction,
        solid_solutions,
        kinetic_species

    export ChemicalState,
        PhaseQuantities,
        temperature,
        pressure,
        set_temperature!,
        set_pressure!,
        set_neutral_pH!,
        moles,
        set_quantity!,
        rescale!,
        mass,
        volume,
        pH,
        pOH,
        porosity,
        saturation,
        volume_fractions,
        chemical_shrinkage,
        missing_molar_volumes

    export extract_primary_species

    export read_thermofun_database,
        build_species,
        build_reactions,
        build_solid_solutions,
        get_compatible_species,
        HKF_SI_CONVERSIONS

    export merge_json

    export AbstractActivityModel,
        DiluteSolutionModel,
        HKFActivityModel,
        DaviesActivityModel,
        activity_model,
        build_potentials,
        REJ_HKF,
        REJ_CHARGE_DEFAULT,
        hkf_debye_huckel_params

    export EquilibriumProblem

    export EquilibriumSolver,
        equilibrate,
        _DEFAULT_SOLVER_FACTORY

    export StateView,
        KineticFunc,
        KINETICS_RATE_MODELS,
        KINETICS_RATE_FACTORIES,
        add_kinetics_rate_model,
        arrhenius_rate_constant,
        saturation_ratio,
        RateModelCatalyst,
        RateMechanism,
        parrot_killoh,
        PK_PARAMS_C3S,
        PK_PARAMS_C2S,
        PK_PARAMS_C3A,
        PK_PARAMS_C4AF,
        parrot_killoh_avrami,
        PK84_PARAMS_C3S,
        PK84_PARAMS_C2S,
        PK84_PARAMS_C3A,
        PK84_PARAMS_C4AF,
        waller,
        WALLER_PARAMS_FLY_ASH,
        WALLER_PARAMS_SILICA_FUME,
        WALLER_PARAMS_SLAG,
        blaine_factor,
        humidity_factor,
        powers_alpha_max

    export AbstractSurfaceModel,
        FixedSurfaceArea,
        BETSurfaceArea,
        surface_area,
        transition_state,
        first_order_rate,
        KineticReaction,
        molar_mass

    export KineticsProblem,
        build_u0,
        build_kinetics_params,
        build_kinetics_ode

    export KineticsSolver,
        integrate,
        _DEFAULT_KINETICS_SOLVER_FACTORY

    export reaction_extents,
        extent_residual,
        state_at,
        speciated_states,
        degrees_of_hydration,
        mean_degree_of_hydration

    export AbstractCalorimeter,
        IsothermalCalorimeter,
        SemiAdiabaticCalorimeter,
        n_extra_states,
        extend_u0,
        extend_ode!,
        heat_rate,
        heat_flow,
        cumulative_heat,
        temperature_profile

    function __init__()
        for (k, v) in THERMO_MODELS
            THERMO_FACTORIES[k] = build_thermo_factories(v)
        end
        for (k, v) in KINETICS_RATE_MODELS
            KINETICS_RATE_FACTORIES[k] = _build_kinetics_rate_factory(v)
        end
        return
    end

end

module ChemistryLab

import Base: ==, +, -, *, /, //
using Crayons
using CSV
using DataFrames
using JSON3
using LinearAlgebra
using OrderedCollections
using PeriodicTable
using PrettyTables
using ProgressBars
using Unicode
import Unitful: u, g, cm, K, J, mol, bar, Quantity, uconvert, ustrip, unit, uparse, upreferred, preferunits, @u_str

preferunits(u"g, cm, K, mol, s"...)

const COL_CHARGE = crayon"cyan bold"
const COL_PAR = crayon"magenta bold"
const COL_STOICH_INT = crayon"red bold"
const COL_STOICH_EXT = crayon"yellow bold"

function print_title(title; crayon=:none, indent="", style=:none)
    draw = crayon !=:none ? x->println(crayon(x)) : println
    if style == :underline
        width = length(title)
        draw(indent * "$title")
        draw(indent * "─"^width)
    elseif style == :box
        width = length(title) + 4
        draw(indent * "┌" * "─"^width * "┐")
        draw(indent * "│  $title  │")
        draw(indent * "└" * "─"^width * "┘")
    else
        draw(indent * "$title") 
    end
    print(crayon"reset")
end

include("parsing_utils.jl")
include("thermodynamical_functions.jl")
include("formulas.jl")
include("species.jl")
include("databases.jl")
include("stoich_matrix.jl")
include("reactions.jl")

export ATOMIC_ORDER, OXIDE_ORDER, stoich_coef_round, phreeqc_to_unicode, colored_formula, unicode_to_phreeqc, parse_formula, extract_charge, calculate_molar_mass, to_mendeleev, parse_equation, format_equation, colored_equation
export Callable, ThermoFunction
export Formula, expr, phreeqc, unicode, colored, composition, charge, apply, check_mendeleev
export AggregateState, AS_UNDEF, AS_AQUEOUS, AS_CRYSTAL, AS_GAS
export Class, SC_UNDEF, SC_AQSOLVENT, SC_AQSOLUTE, SC_COMPONENT, SC_GAS_FLUID
export AbstractSpecies, Species, CemSpecies
export name, symbol, formula, mendeleev_filter, cemformula, atoms, atoms_charge, oxides, oxides_charge, components, aggregate_state, class, properties
export merge_json, read_thermofun, extract_primary_species
export union_atoms, print_stoich_matrix, canonical_stoich_matrix, stoich_matrix, stoich_matrix_to_equations, stoich_matrix_to_reactions
export Reaction, CemReaction, reactants, products, simplify_reaction
@eval export $(Symbol.(EQUAL_OPS)...) 

end

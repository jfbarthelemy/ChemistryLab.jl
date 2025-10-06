module CementChemistry

import Base: ==, +, -, *, /, //
using LinearAlgebra
using OrderedCollections
using JSON3, CSV, DataFrames
using Unicode
import Unitful: u, g, cm, K, J, mol, Quantity, uconvert
using PeriodicTable
using PhysicalConstants.CODATA2022
using PrettyTables, Crayons

const COL_CHARGE = crayon"cyan bold"
const COL_PAR = crayon"magenta bold"
const COL_STOICH_INT = crayon"red bold"
const COL_STOICH_EXT = crayon"yellow bold"

include("parsing_utils.jl")
include("formulas.jl")
include("species.jl")
include("databases.jl")
include("stoich_matrix.jl")
include("reactions.jl")

export ATOMIC_ORDER, OXIDE_ORDER, stoich_coef_round, phreeqc_to_unicode, colored_formula, unicode_to_phreeqc, parse_formula, extract_charge, calculate_molar_mass, to_mendeleev, parse_equation, format_equation, colored_equation
export Formula, expr, phreeqc, unicode, colored, composition, charge, apply
export AbstractSpecies, Species, CemSpecies
export name, symbol, formula, cemformula, atoms, atoms_charge, oxides, oxides_charge, components, properties
export merge_json, parse_cemdata18_thermofun, extract_primary_species
export union_atoms, print_stoich_matrix, canonical_stoich_matrix, stoich_matrix, stoich_matrix_to_equations, stoich_matrix_to_reactions
export Reaction, CemReaction, reactants, products
@eval export $(Symbol.(EQUAL_OPS)...) 

end

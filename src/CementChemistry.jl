module CementChemistry

export ATOMIC_ORDER, OXIDE_ORDER, stoich_coef_round, phreeqc_to_unicode, unicode_to_phreeqc, parse_formula, extract_charge, calculate_molar_mass, to_mendeleev, parse_equation, format_equation
export Formula, expr, phreeqc, unicode, composition, charge
export AbstractSpecies, Species, CemSpecies
export name, symbol, formula, cemformula, atoms, atoms_charge, oxides, oxides_charge, components, properties
export merge_json, parse_cemdata18_thermofun, extract_primary_species
export union_atoms, stoich_matrix, stoich_matrix_to_equations

import Base: ==, +, -, *, /, //
using LinearAlgebra
using JSON, JSON3, CSV, DataFrames
using Unicode
import Unitful: u, g, cm, K, J, mol, Quantity, uconvert
using PeriodicTable
using PhysicalConstants.CODATA2022
using PrettyTables

include("parsing_utils.jl")
include("formulas.jl")
include("species.jl")
include("databases.jl")
include("stoich_matrix.jl")




end

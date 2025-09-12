module CementChemistry

export phreeqc_to_unicode, unicode_to_phreeqc, parse_formula, extract_charge, calculate_molar_mass
export Species, name, formula, symbol, phreeqc, unicode, atoms, charge, molar_mass, properties

import Base: ==, +, -, *, /, //
using Unicode
import Unitful: u, g, cm, K, J, mol, Quantity, uconvert
using PeriodicTable
using PhysicalConstants.CODATA2022

include("parsing_utils.jl")
include("species.jl")




end

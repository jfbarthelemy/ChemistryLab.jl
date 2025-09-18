# CementChemistry [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jfbarthelemy.github.io/CementChemistry.jl/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jfbarthelemy.github.io/CementChemistry.jl/dev/) [![Build Status](https://github.com/jfbarthelemy/CementChemistry.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jfbarthelemy/CementChemistry.jl/actions/workflows/CI.yml?query=branch%3Amain)

CementChemistry.jl is a Julia library for chemical modeling of cements and aqueous solutions, supporting ThermoFun and Cemdata databases.

## Features

- **Chemical formula handling**: Create, convert, and display formulas with charge management and Unicode/Phreeqc notation.
- **Chemical species management**: `Species` and `CemSpecies` types to represent solution and solid phase species.
- **Stoichiometric matrices**: Automatic construction of matrices for reaction and equilibrium analysis.
- **Database interoperability**: Import and merge ThermoFun (.json) and Cemdata (.dat) data.
- **Parsing tools**: Convert chemical notations, extract charges, calculate molar mass, and more.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/jfbarthelemy/CementChemistry.jl")
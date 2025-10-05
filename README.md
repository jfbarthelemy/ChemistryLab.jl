# CementChemistry [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jfbarthelemy.github.io/CementChemistry.jl/dev/) [![Build Status](https://github.com/jfbarthelemy/CementChemistry.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jfbarthelemy/CementChemistry.jl/actions/workflows/CI.yml?query=branch%3Amain)
<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jfbarthelemy.github.io/CementChemistry.jl/stable/)  -->

CementChemistry.jl is a computational chemistry toolkit for modeling low-carbon cementitious materials and aqueous solutions. Designed for researchers, engineers, and developers working with cement chemistry, it provides formula handling, species management, stoichiometric matrix construction, and database interoperability (ThermoFun and Cemdata). Main features include chemical formula parsing, Unicode/Phreeqc notation conversion, reaction and equilibrium analysis, and data import/export.

<!-- ## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Prerequisites](#prerequisites)
- [Contributing](#contributing)
- [Known Bugs](#known-bugs)
- [FAQ](#faq)
- [License](#license)
- [Credits and Acknowledgements](#credits-and-acknowledgements)
- [Contact](#contact) -->

## Features

- **Chemical formula handling**: Create, convert, and display formulas with charge management and Unicode/Phreeqc notation.
- **Chemical species management**: `Species` and `CemSpecies` types to represent solution and solid phase species.
- **Stoichiometric matrices**: Automatic construction of matrices for reaction and equilibrium analysis.
- **Database interoperability**: Import and merge ThermoFun (.json) and Cemdata (.dat) data.
- **Parsing tools**: Convert chemical notations, extract charges, calculate molar mass, and more.

## Prerequisites

- Julia 1.6 or higher
- Dependencies: JSON, JSON3, CSV, DataFrames, Unicode, Unitful, PeriodicTable, PhysicalConstants.CODATA2022, PrettyTables
- No API account or special access required 

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/jfbarthelemy/CementChemistry.jl")
```

## Project Status

- Version: x.y.z
- Last update: 2025-09-25
- Status: In development / Stable

## Usage

See the [documentation and tutorials](https://jfbarthelemy.github.io/CementChemistry.jl/dev/
) for examples on formula creation, species management, reaction parsing, and database merging.

## Contributing

Contributions are welcome! Please submit issues or pull requests via GitHub.

## Known Bugs
See GitHub Issues for current bug reports.

## FAQ

Frequently asked questions are addressed in the [documentation](https://jfbarthelemy.github.io/CementChemistry.jl/dev/).

## License

MIT License. See [LICENSE](LICENSE) for details.

## Credits and Acknowledgements

Developed by Jean-François Barthélémy and Anthony Soive. Thanks to contributors and the open-source community.

## Contact

For questions or support, contact the authors via GitHub or email listed in the repository.
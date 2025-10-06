# Database Interoperability

So far, we have looked at the possibility of creating and manipulating any species, whether they exist or not. Creating an H₂O⁺⁴ molecule, for example, is not a problem.

```@setup database_interoperability
    using CementChemistry
```

```@example database_interoperability
HSO4 = Species("H₂O⁴⁺")
```

However, you will admit that it is a little strange...

## Cemdata18 and PSI-Nagra-12-07 Databases

This is why Cement Chemistry relies on existing databases, in particular [Cemdata18](https://www.empa.ch/web/s308/thermodynamic-data) and [PSI-Nagra-12-07](https://www.psi.ch/en/les/thermodynamic-databases). Cemdata18 is a chemical thermodynamic database for hydrated Portland cements and alkali-activated materials. PSI-Nagra is a Chemical Thermodynamic Database. The formalism adopted for these databases is that of [Thermofun](https://thermohub.org/thermofun/thermofun/) which is a universal open-source client that delivers thermodynamic properties of substances and reactions at the temperature and pressure of interest. The information is stored in json files.

With Cementchemistry, you can parse a ThermoFun-like json file and return DataFrames for:

- elements:

```@example database_interoperability
df_elements, df_substances, df_reactions = parse_cemdata18_thermofun("../../../data/cemdata18-merged.json")
show(df_elements, allcols=true, allrows=true)
```
[`CementChemistry.parse_cemdata18_thermofun`](@ref)

- species (aqueous, solid or gaseous phases):

```@example database_interoperability
show(df_substances, allcols=true, allrows=false)
```

- reactions contained in the database
```@example database_interoperability
show(df_reactions, allcols=true, allrows=false)
```

## Primary species extraction

It is also possible to retrieve primary species from the Cemdata18 database, primary species being the designation of a subset of species for which any species can be represented as the linear combination of primary species.

```@example database_interoperability
df_primaries = extract_primary_species("../../../data/CEMDATA18-31-03-2022-phaseVol.dat")
show(df_primaries, allcols=true, allrows=true)
```
[`CementChemistry.extract_primary_species`](@ref)

---


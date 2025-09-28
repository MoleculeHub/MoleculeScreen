# MoleculeScreen.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](#installation)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/JuliaDiff/BlueStyle)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

A Julia package for molecular filtering and screening using both property-based and SMARTS-based filters. 

## Features

- **Property-based filters**: Lipinski's Rule of Five, Veber's rules, Ghose filter, Egan filter, Muegge filter, Pfizer's 3/75 rule, GSK's 4/400 rule, Golden Triangle
- **SMARTS-based filters**: PAINS, Brenk, BMS, SureChEMBL, MLSMR, Inpharmatica, LINT, Eli Lilly, and Glaxo filters
- **High-level screening functions**: Screen individual molecules or filter collections

## Installation

```julia
using Pkg
Pkg.add("MoleculeScreen")  
```

## Quick Start

```julia
using MoleculeScreen
using MoleculeFlow

# Create a molecule
mol = mol_from_smiles("CCO")  

# Apply property filters
results = apply_property_filters(mol; filters=[:lipinski, :veber])
println(results)  # Dict(:lipinski => true, :veber => true)

# Apply SMARTS filters
smarts_results = apply_smarts_filters(mol; filters=["pains", "brenk"])
println(smarts_results)  # Dict("pains" => true, "brenk" => true)

# Comprehensive screening
screening_results = screen_molecule(mol)
println(screening_results["overall_pass"])  # true
```

## Property-Based Filters

### Available Filters

- `:lipinski` - Lipinski's Rule of Five (MW ≤ 500, LogP ≤ 5, HBD ≤ 5, HBA ≤ 10)
- `:veber` - Veber's rules (Rotatable bonds ≤ 10, TPSA ≤ 140)
- `:ghose` - Ghose filter (MW: 160-480, LogP: -0.4 to 5.6, etc.)
- `:egan` - Egan's filter (LogP ≤ 5.88, TPSA ≤ 131)
- `:muegge` - Muegge filter (MW: 200-600, LogP: -2 to 5, etc.)
- `:pfizer` - Pfizer's 3/75 rule (flags LogP > 3 AND TPSA < 75)
- `:gsk` - GSK's 4/400 rule (HBD ≤ 4, MW ≤ 400)
- `:golden_triangle` - Golden Triangle (LogP: -0.5 to 4.5, MW: 200-500)

### Usage

```julia
# Single filter
is_drug_like = lipinski_ro5(mol)

# Multiple filters
results = apply_property_filters(mol; filters=[:lipinski, :veber, :egan])
passes_all = all(values(results))
```

## SMARTS-Based Filters

### Available Filter Sets

The package loads SMARTS patterns from `data/rules_data.csv` with the following major filter sets:

- `"pains"` 
- `"elililly"`
- `"surechembl"`
- `"mlsmr"`
- `"bms"`
- `"brenk"` 
- `"inpharmatica"` 
- `"glaxo"` 
- `"lint"` 

### Usage

```julia
available = get_available_smarts_filters()

# Single filter check
is_clean = check_smarts_filter(mol, "pains")

# Multiple filters
results = apply_smarts_filters(mol; filters=["pains", "brenk"])

# Get specific violations
violations = get_smarts_violations(mol, "pains")

is_clean = is_clean_molecule(mol; smarts_filters=["pains", "brenk"])
```

## High-Level Screening

### Screen Individual Molecules

```julia
# Comprehensive screening with default filters
results = screen_molecule(mol)

# Custom screening
results = screen_molecule(mol;
    property_filters=[:lipinski, :veber],
    smarts_filters=["pains", "brenk"])
```

### Filter Molecule Collections

```julia
molecules = [
    mol_from_smiles("CCO"),           
    mol_from_smiles("c1ccccc1N"),     
    mol_from_smiles("CC(=O)O"),       
    missing                           
]

# Filter collection (returns only passing molecules)
filtered = filter_molecules(molecules;
    property_filters=[:lipinski],
    smarts_filters=["pains"])

# Get detailed results as DataFrame
detailed_results = filter_molecules(molecules; return_detailed=true)
```

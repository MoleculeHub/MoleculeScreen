module MoleculeScreen

using MoleculeFlow
using CSV
using DataFrames

include("property_rules.jl")
include("smarts_rules.jl")

export lipinski_ro5, veber_rules, ghose_filter, egan_filter, muegge_filter
export pfizer_3_75_rule, gsk_4_400_rule, golden_triangle, apply_property_filters

export check_smarts_filter, apply_smarts_filters, get_smarts_violations
export is_clean_molecule, get_available_smarts_filters, get_smarts_violations_with_descriptions

"""
    screen_molecules(mol::Molecule;
                    property_filters::Vector{Symbol}=[:lipinski],
                    smarts_filters::Vector{String}=["pains"]) -> Dict

Comprehensive screening of a molecule using both property-based and SMARTS-based filters.

# Arguments

  - `mol::Molecule`: Molecule to screen
  - `property_filters::Vector{Symbol}`: Property-based filters to apply
  - `smarts_filters::Vector{String}`: SMARTS-based filters to apply

# Returns

  - `Dict`: Results containing both property and SMARTS filter results, plus overall pass/fail

# Examples

```julia
mol = mol_from_smiles("CCO")  # Ethanol
results = screen_molecules(mol)
# Returns comprehensive screening results

# Custom filtering
results = screen_molecules(
    mol; property_filters = [:lipinski, :veber], smarts_filters = ["pains", "brenk"]
)
```
"""
function screen_molecules(
    mol::Molecule;
    property_filters::Vector{Symbol} = [:lipinski],
    smarts_filters::Vector{String} = ["pains"],
)
    property_results = apply_property_filters(mol; filters = property_filters)

    smarts_results = apply_smarts_filters(mol; filters = smarts_filters)

    property_pass = all(values(property_results))
    smarts_pass = all(r["pass"] for r in values(smarts_results))
    overall_pass = property_pass && smarts_pass

    return Dict(
        "property_filters" => property_results,
        "smarts_filters" => smarts_results,
        "property_pass" => property_pass,
        "smarts_pass" => smarts_pass,
        "overall_pass" => overall_pass,
        "molecule_valid" => mol.valid,
    )
end

"""
    screen_molecules(mols::Vector{Union{Molecule,Missing}};
                    property_filters::Vector{Symbol}=[:lipinski],
                    smarts_filters::Vector{String}=["pains"]) -> Vector{Dict}

Comprehensive screening of multiple molecules using both property-based and SMARTS-based filters.

# Arguments

  - `mols::Vector{Union{Molecule,Missing}}`: Vector of molecules to screen
  - `property_filters::Vector{Symbol}`: Property-based filters to apply
  - `smarts_filters::Vector{String}`: SMARTS-based filters to apply

# Returns

  - `Vector{Dict}`: Vector of screening results for each molecule

# Examples

```julia
mols = [mol_from_smiles("CCO"), mol_from_smiles("c1ccccc1N")]
results = screen_molecules(mols; smarts_filters = ["pains"])
```
"""
function screen_molecules(
    mols::Vector{T} where {T <: Union{Molecule, Missing}};
    property_filters::Vector{Symbol} = [:lipinski],
    smarts_filters::Vector{String} = ["pains"],
)
    results = []

    for (i, mol) in enumerate(mols)
        if ismissing(mol)
            push!(results, Dict(
                "index" => i,
                "molecule_valid" => missing,
                "property_pass" => missing,
                "smarts_pass" => missing,
                "overall_pass" => missing,
                "property_filters" => missing,
                "smarts_filters" => missing
            ))
            continue
        end

        result = screen_molecules(
            mol; property_filters = property_filters, smarts_filters = smarts_filters
        )
        result["index"] = i
        push!(results, result)
    end

    return results
end

export screen_molecules

end

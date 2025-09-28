module MoleculeScreen

using MoleculeFlow
using CSV
using DataFrames

include("property_rules.jl")
include("smarts_rules.jl")

export lipinski_ro5, veber_rules, ghose_filter, egan_filter, muegge_filter
export pfizer_3_75_rule, gsk_4_400_rule, golden_triangle, apply_property_filters

export check_smarts_filter, apply_smarts_filters, get_smarts_violations
export is_clean_molecule, get_available_smarts_filters

"""
    screen_molecule(mol::Molecule;
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
results = screen_molecule(mol)
# Returns comprehensive screening results

# Custom filtering
results = screen_molecule(mol;
    property_filters=[:lipinski, :veber],
    smarts_filters=["pains", "brenk"])
```
"""
function screen_molecule(mol::Molecule;
                        property_filters::Vector{Symbol}=[:lipinski],
                        smarts_filters::Vector{String}=["pains"])

    property_results = apply_property_filters(mol; filters=property_filters)

    smarts_results = apply_smarts_filters(mol; filters=smarts_filters)

    property_pass = all(values(property_results))
    smarts_pass = all(values(smarts_results))
    overall_pass = property_pass && smarts_pass

    return Dict(
        "property_filters" => property_results,
        "smarts_filters" => smarts_results,
        "property_pass" => property_pass,
        "smarts_pass" => smarts_pass,
        "overall_pass" => overall_pass,
        "molecule_valid" => mol.valid
    )
end

"""
    filter_molecules(mols::Vector{Union{Molecule,Missing}};
                    property_filters::Vector{Symbol}=[:lipinski],
                    smarts_filters::Vector{String}=["pains"],
                    return_detailed::Bool=false) -> Union{Vector{Union{Molecule,Missing}}, DataFrame}

Filter a vector of molecules using specified criteria.

# Arguments
- `mols::Vector{Union{Molecule,Missing}}`: Vector of molecules to filter
- `property_filters::Vector{Symbol}`: Property-based filters to apply
- `smarts_filters::Vector{String}`: SMARTS-based filters to apply
- `return_detailed::Bool`: If true, return DataFrame with detailed results; if false, return filtered molecules

# Returns
- `Vector{Union{Molecule,Missing}}`: Filtered molecules (if return_detailed=false)
- `DataFrame`: Detailed filtering results (if return_detailed=true)

# Examples
```julia
mols = [mol_from_smiles("CCO"), mol_from_smiles("c1ccccc1N")]
filtered = filter_molecules(mols; smarts_filters=["pains"])

# Get detailed results
results_df = filter_molecules(mols; return_detailed=true)
```
"""
function filter_molecules(mols::Vector{T} where T<:Union{Molecule,Missing};
                         property_filters::Vector{Symbol}=[:lipinski],
                         smarts_filters::Vector{String}=["pains"],
                         return_detailed::Bool=false)

    results = []

    for (i, mol) in enumerate(mols)
        if ismissing(mol)
            push!(results, (index=i, mol=mol, result=missing))
            continue
        end

        result = screen_molecule(mol;
                               property_filters=property_filters,
                               smarts_filters=smarts_filters)
        push!(results, (index=i, mol=mol, result=result))
    end

    if return_detailed
        return DataFrame([
            (index=r.index,
             molecule_valid=ismissing(r.result) ? missing : r.result["molecule_valid"],
             property_pass=ismissing(r.result) ? missing : r.result["property_pass"],
             smarts_pass=ismissing(r.result) ? missing : r.result["smarts_pass"],
             overall_pass=ismissing(r.result) ? missing : r.result["overall_pass"])
            for r in results
        ])
    else
        passing_mols = Union{Molecule,Missing}[]
        for r in results
            if ismissing(r.result) || r.result["overall_pass"]
                push!(passing_mols, r.mol)
            end
        end
        return passing_mols
    end
end

export screen_molecule, filter_molecules

end
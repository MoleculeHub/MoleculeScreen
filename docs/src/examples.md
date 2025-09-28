# Examples

## Drug Screening Examples

```julia
using MoleculeScreen
using MoleculeFlow

# Test well-known drugs
drugs = [
    ("CCO", "Ethanol"),
    ("CC(C)C1=CC(=C(C=C1)C(C)C)C(=O)C(=O)O", "Ibuprofen"),
    ("CC(=O)OC1=CC=CC=C1C(=O)O", "Aspirin"),
    ("C=CC(=O)N1CCC[C@H](C1)N2C3=C(C(=N2)C4=CC=C(C=C4)OC5=CC=CC=C5)C(=NC=N3)N", "Ibrutinib")
]

for (smiles, name) in drugs
    mol = mol_from_smiles(smiles)

    if !mol.valid
        println("$name: Invalid SMILES")
        continue
    end

    results = screen_molecules(mol;
        property_filters=[:lipinski, :veber],
        smarts_filters=["pains", "brenk"])

    println("\n$name:")
    println("Overall pass:", results["overall_pass"])
    println("Property pass:", results["property_pass"])
    println("SMARTS pass:", results["smarts_pass"])

    # Show any violations
    for (filter_name, filter_result) in results["smarts_filters"]
        if !filter_result["pass"]
            println("$filter_name violations:")
            for violation in filter_result["violations"]
                println("- $(violation.description)")
            end
        end
    end
end
```

## Error Handling Examples

### Robust Pipeline

```julia
function production_screening(smiles_list; log_errors=true)
    results = []
    errors = []

    for (i, smiles) in enumerate(smiles_list)
        try
            mol = mol_from_smiles(smiles)

            if !mol.valid
                error_msg = "Invalid SMILES: $smiles"
                if log_errors
                    push!(errors, (index=i, smiles=smiles, error=error_msg))
                end
                continue
            end

            result = screen_molecules(mol;
                property_filters=[:lipinski, :veber],
                smarts_filters=["pains", "brenk"])

            result["index"] = i
            result["smiles"] = smiles
            push!(results, result)

        catch e
            error_msg = "Processing error for $smiles: $(string(e))"
            if log_errors
                push!(errors, (index=i, smiles=smiles, error=error_msg))
            end
            println("Warning: $error_msg")
        end
    end

    return (results=results, errors=errors)
end

test_smiles = ["CCO", "INVALID", "c1ccccc1", "ERROR_PRONE_SMILES"]
screening_output = production_screening(test_smiles)

println("Successfully screened: $(length(screening_output.results)) molecules")
println("Errors encountered: $(length(screening_output.errors))")

for error in screening_output.errors
    println("Error at index $(error.index): $(error.error)")
end
```
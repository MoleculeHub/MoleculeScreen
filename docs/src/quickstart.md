# Quickstart

```julia
using MoleculeScreen
using MoleculeFlow

# Create a molecule from SMILES
mol = mol_from_smiles("CCO")  

!mol.valid

# Screen with default filters
results = screen_molecules(mol)
println("Molecule passes screening: ", results["overall_pass"])
```

# Basic Screening 
```julia
# Single Property for a single Molecule 
mol = mol_from_smiles("CC(C)C1=CC(=C(C=C1)C(C)C)C(=O)C(=O)O")  # Ibuprofen

if mol.valid
    # Test individual filters
    println("Lipinski: ", lipinski_ro5(mol))
    println("Veber: ", veber_rules(mol))
    println("Ghose: ", ghose_filter(mol))
end
```

# Advanced Screening 
```julia
results = screen_molecules(mol;
    property_filters=[:lipinski, :veber],
    smarts_filters=["pains", "brenk"])

# Result structure:
# {
#   "property_filters" => {:lipinski => true, :veber => true},
#   "smarts_filters" => {
#     "pains" => {"pass" => true, "violations" => []},
#     "brenk" => {"pass" => false, "violations" => [(smarts="...", description="...")]}
#   },
#   "property_pass" => true,
#   "smarts_pass" => false,
#   "overall_pass" => false,
#   "molecule_valid" => true
# }
```
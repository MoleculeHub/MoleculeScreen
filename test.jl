using MoleculeScreen
using MoleculeFlow

mol = mol_from_smiles("CC")
results = screen_molecule(mol)

# Custom screening
results = screen_molecule(
    mol; property_filters = [:lipinski, :veber], smarts_filters = ["pains", "brenk"]
)

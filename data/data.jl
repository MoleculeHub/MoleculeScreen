using Pkg
Pkg.add("TidierData")
using TidierData
using CSV
using MoleculeFlow

df = DataFrame(CSV.File("./data/rules_data.csv"))

mol = mol_from_smiles(
    "C=CC(=O)N1CCC[C@H](C1)N2C3=C(C(=N2)C4=CC=C(C=C4)OC5=CC=CC=C5)C(=NC=N3)N"
)

for (pat, desc) in zip(df_eli.smarts, df_eli.description)
    try
        if has_substructure_match(mol, pat) == true
            println("$pat $desc")
            # highlight_substructure(mol, pat)
        end
    catch e
        # println("$pat $desc")
    end
end

highlight_substructure(mol, "O-c1ccccc1")

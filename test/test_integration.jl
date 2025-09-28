using Test
using MoleculeScreen
using MoleculeFlow

@testset "Integration tests" begin
    ethanol = mol_from_smiles("CCO")
    aspirin = mol_from_smiles("CC(=O)OC1=CC=CC=C1C(=O)O")
    aniline = mol_from_smiles("c1ccccc1N")
    large_mol = mol_from_smiles("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC")

    @testset "Full molecule screening" begin
        results = screen_molecule(ethanol)
        @test results["molecule_valid"] == true
        @test results["property_pass"] == true
        @test results["smarts_pass"] == true
        @test results["overall_pass"] == true

        results_aniline = screen_molecule(aniline)
        @test results_aniline["molecule_valid"] == true
        @test results_aniline["property_pass"] == true
        @test results_aniline["smarts_pass"] isa Bool
        @test results_aniline["overall_pass"] isa Bool

        results_custom = screen_molecule(
            aspirin;
            property_filters = [:lipinski, :veber],
            smarts_filters = ["pains", "brenk"],
        )
        @test haskey(results_custom["property_filters"], :lipinski)
        @test haskey(results_custom["property_filters"], :veber)
        @test haskey(results_custom["smarts_filters"], "pains")
        @test haskey(results_custom["smarts_filters"], "brenk")
    end

    @testset "Bulk molecule filtering" begin
        molecules = Union{Molecule, Missing}[ethanol, aspirin, aniline, large_mol]

        filtered = filter_molecules(molecules)
        @test length(filtered) <= length(molecules)

        detailed = filter_molecules(molecules; return_detailed = true)
        @test nrow(detailed) == length(molecules)
        @test hasproperty(detailed, :index)
        @test hasproperty(detailed, :molecule_valid)
        @test hasproperty(detailed, :property_pass)
        @test hasproperty(detailed, :smarts_pass)
        @test hasproperty(detailed, :overall_pass)

        filtered_custom = filter_molecules(
            molecules; property_filters = [:lipinski], smarts_filters = ["pains"]
        )
        @test length(filtered_custom) >= 1
    end

    @testset "Handle missing molecules" begin
        molecules_with_missing = Union{Molecule, Missing}[ethanol, missing, aspirin]

        filtered = filter_molecules(molecules_with_missing)
        @test any(ismissing, filtered)

        detailed = filter_molecules(molecules_with_missing; return_detailed = true)
        @test nrow(detailed) == 3
        @test ismissing(detailed[2, :molecule_valid])
    end
end

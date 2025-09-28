using Test
using MoleculeScreen
using MoleculeFlow

@testset "Integration tests" begin
    ethanol = mol_from_smiles("CCO")
    aspirin = mol_from_smiles("CC(=O)OC1=CC=CC=C1C(=O)O")
    aniline = mol_from_smiles("c1ccccc1N")
    large_mol = mol_from_smiles("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC")

    @testset "Full molecule screening" begin
        results = screen_molecules(ethanol)
        @test results["molecule_valid"] == true
        @test results["property_pass"] == true
        @test results["smarts_pass"] == true
        @test results["overall_pass"] == true

        results_aniline = screen_molecules(aniline)
        @test results_aniline["molecule_valid"] == true
        @test results_aniline["property_pass"] == true
        @test results_aniline["smarts_pass"] isa Bool
        @test results_aniline["overall_pass"] isa Bool

        results_custom = screen_molecules(
            aspirin;
            property_filters = [:lipinski, :veber],
            smarts_filters = ["pains", "brenk"],
        )
        @test haskey(results_custom["property_filters"], :lipinski)
        @test haskey(results_custom["property_filters"], :veber)
        @test haskey(results_custom["smarts_filters"], "pains")
        @test haskey(results_custom["smarts_filters"], "brenk")
        @test haskey(results_custom["smarts_filters"]["pains"], "pass")
        @test haskey(results_custom["smarts_filters"]["pains"], "violations")
    end

    @testset "Bulk molecule screening" begin
        molecules = Union{Molecule, Missing}[ethanol, aspirin, aniline, large_mol]

        results = screen_molecules(molecules)
        @test length(results) == length(molecules)

        for result in results
            @test haskey(result, "index")
            @test haskey(result, "molecule_valid")
            @test haskey(result, "property_pass")
            @test haskey(result, "smarts_pass")
            @test haskey(result, "overall_pass")
        end

        results_custom = screen_molecules(
            molecules; property_filters = [:lipinski], smarts_filters = ["pains"]
        )
        @test length(results_custom) == length(molecules)
    end

    @testset "Handle missing molecules" begin
        molecules_with_missing = Union{Molecule, Missing}[ethanol, missing, aspirin]

        results = screen_molecules(molecules_with_missing)
        @test length(results) == 3
        @test ismissing(results[2]["molecule_valid"])
        @test results[2]["index"] == 2
    end
end

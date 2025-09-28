using Test
using MoleculeScreen
using MoleculeFlow

@testset "Data loading tests" begin
    @testset "SMARTS rules loading" begin
        available_filters = get_available_smarts_filters()
        @test length(available_filters) > 0

        @test "pains" in available_filters
        @test "brenk" in available_filters
        @test "elililly" in available_filters
        @test "surechembl" in available_filters

        ethanol = mol_from_smiles("CCO")
        @test check_smarts_filter(ethanol, "pains") isa Bool

        @test_logs (:warn, r"Unknown SMARTS filter") check_smarts_filter(
            ethanol, "nonexistent_filter"
        )
    end

    @testset "Rule set sizes" begin
        available = get_available_smarts_filters()

        if "pains" in available
            pains_violations = get_smarts_violations(mol_from_smiles("c1ccccc1N"), "pains")
            @test pains_violations isa Vector{String}
        end
    end
end

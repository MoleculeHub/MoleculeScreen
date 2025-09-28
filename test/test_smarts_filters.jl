using Test
using MoleculeScreen
using MoleculeFlow

@testset "SMARTS-based filters" begin
    ethanol = mol_from_smiles("CCO")
    aniline = mol_from_smiles("c1ccccc1N")  
    catechol = mol_from_smiles("c1ccc(O)c(O)c1")  

    @testset "Available filters" begin
        available = get_available_smarts_filters()
        @test "pains" in available
        @test "brenk" in available
        @test length(available) > 5
    end

    @testset "Individual SMARTS checks" begin
        @test check_smarts_filter(ethanol, "pains") == true  
        pains_result = check_smarts_filter(aniline, "pains")
        @test pains_result isa Bool  
        catechol_result = check_smarts_filter(catechol, "pains")
        @test catechol_result isa Bool  
    end

    @testset "Case insensitive filter names" begin
        @test check_smarts_filter(ethanol, "PAINS") == true
        @test check_smarts_filter(ethanol, "Pains") == true
        @test check_smarts_filter(ethanol, "pains") == true
    end

    @testset "SMARTS violations" begin
        violations = get_smarts_violations(aniline, "pains")
        @test violations isa Vector{String}  

        violations_clean = get_smarts_violations(ethanol, "pains")
        @test violations_clean isa Vector{String}  
        @test length(violations_clean) == 0  
    end

    @testset "Multiple SMARTS filters" begin
        filters = ["pains", "brenk"]
        results = apply_smarts_filters(ethanol; filters=filters)

        @test haskey(results, "pains")
        @test haskey(results, "brenk")
        @test all(values(results))  

        results_aniline = apply_smarts_filters(aniline; filters=filters)
        @test haskey(results_aniline, "pains")  
    end

    @testset "Clean molecule check" begin
        @test is_clean_molecule(ethanol; smarts_filters=["pains"]) == true
        aniline_result = is_clean_molecule(aniline; smarts_filters=["pains"])
        @test aniline_result isa Bool
    end
end
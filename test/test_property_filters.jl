using Test
using MoleculeScreen
using MoleculeFlow

@testset "Property-based filters" begin
    ethanol = mol_from_smiles("CCO")
    aspirin = mol_from_smiles("CC(=O)OC1=CC=CC=C1C(=O)O")
    caffeine = mol_from_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")

    @testset "Lipinski Rule of Five" begin
        @test lipinski_ro5(ethanol) == true
        @test lipinski_ro5(aspirin) == true
        @test lipinski_ro5(caffeine) == true
    end

    @testset "Veber Rules" begin
        @test veber_rules(ethanol) == true
        @test veber_rules(aspirin) == true
        @test veber_rules(caffeine) == true
    end

    @testset "Ghose Filter" begin
        @test ghose_filter(ethanol) == false  
        @test ghose_filter(aspirin) isa Bool  
        @test ghose_filter(caffeine) isa Bool  
    end

    @testset "Egan Filter" begin
        @test egan_filter(ethanol) == true
        @test egan_filter(aspirin) == true
        @test egan_filter(caffeine) == true
    end

    @testset "Golden Triangle" begin
        @test golden_triangle(ethanol) == false 
        @test golden_triangle(aspirin) isa Bool 
        @test golden_triangle(caffeine) isa Bool
    end

    @testset "Apply multiple filters" begin
        filters = [:lipinski, :veber, :egan]
        results = apply_property_filters(aspirin; filters=filters)

        @test haskey(results, :lipinski)
        @test haskey(results, :veber)
        @test haskey(results, :egan)
        @test all(values(results))
    end
end
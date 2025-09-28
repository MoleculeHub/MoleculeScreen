using Test
using MoleculeScreen
using MoleculeFlow
using DataFrames
using CSV
using Aqua

@testset "MoleculeScreen.jl" begin
    @testset "Aqua.jl quality assurance" begin
        Aqua.test_all(MoleculeScreen)
    end

    include("test_data_loading.jl")
    include("test_property_filters.jl")
    include("test_smarts_filters.jl")
    include("test_integration.jl")
end

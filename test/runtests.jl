using Test
using MoleculeScreen
using MoleculeFlow
using DataFrames
using CSV

@testset "MoleculeScreen.jl" begin
    include("test_data_loading.jl")
    include("test_property_filters.jl")
    include("test_smarts_filters.jl")
    include("test_integration.jl")
end
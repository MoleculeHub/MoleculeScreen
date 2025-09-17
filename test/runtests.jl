using Test
using MoleculeScreen
using MoleculeFlow
using DataFrames
using CSV

@testset "MoleculeScreen.jl" begin
    include("test_alert_loading.jl")
    include("test_screening.jl")
    include("test_filtering.jl")
    include("test_integration.jl")
end
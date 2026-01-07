using Test
using FSMOSHD
using Serialization

projdir = dirname(dirname(@__FILE__))

@testset "Soil Energy Balance" begin
    include("test_soil_energy_balance.jl")
end

@testset "Regression Tests" begin
    include("test_regression.jl")
end

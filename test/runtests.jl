using Test
using FlexibleSnowModelOSHD
include("testhelpers.jl")

using Documenter

@testset "Doctests" begin
    # Doctest output can drift across Julia versions (type/float printing), and the
    # CI matrix includes the 1.10 compat floor. Pin doctests to the version we
    # develop against; the floor still runs every other test.
    if VERSION >= v"1.11"
        DocMeta.setdocmeta!(
            FlexibleSnowModelOSHD,
            :DocTestSetup,
            :(using FlexibleSnowModelOSHD);
            recursive = true,
        )
        doctest(FlexibleSnowModelOSHD)
    end
end

@testset "Mass Balance" begin
    include("test_mass_balance.jl")
end

@testset "MET Immutability" begin
    include("test_met_immutability.jl")
end

@testset "Soil Energy Balance" begin
    include("test_soil_energy_balance.jl")
end

@testset "Architectures" begin
    include("test_architectures.jl")
end

@testset "SnowSlide" begin
    include("test_snowslide.jl")
end

@testset "SnowTran3D" begin
    include("test_snowtran3d.jl")
end

@testset "Transport step" begin
    include("test_transport_step.jl")
end

@testset "Regression Tests" begin
    include("test_regression.jl")
end

@testset "Mixed land cover" begin
    include("test_mixed_landcover.jl")
end

# GPU test — runs a full CPU-vs-GPU comparison only if CUDA is functional in the (stacked)
# default environment; otherwise it logs a skip and adds no tests. See test_gpu.jl.
include("test_gpu.jl")

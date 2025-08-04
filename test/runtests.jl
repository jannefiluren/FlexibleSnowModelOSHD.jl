using Test
using FSMOSHD

# Absolute path to project folder
projdir = dirname(dirname(@__FILE__))

# Include regression test modules
include("test_data_generator.jl")
include("generate_reference_results.jl") 
include("regression_tests.jl")

# Test data paths
test_data_path = joinpath(projdir, "test", "test_data")
reference_path = joinpath(test_data_path, "reference")

@testset "FSMOSHD Regression Tests" begin

    # Create test data if it doesn't exist
    if !isdir(reference_path) || !isfile(joinpath(reference_path, "reference_results.mat"))
        println("Creating test dataset and reference results...")
        create_test_dataset()
        create_full_reference_dataset()
    end
    
    # Run regression test
    println("Running regression tests...")
    test_passed, max_differences = run_regression_test(test_data_path, tolerance=1e-10)
    
    @test test_passed
    
end
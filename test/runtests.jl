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

    # Configuration matrix: (tile, snfrac)
    configurations = [
        ("open", 0),
        ("open", 3), 
        ("open", 4),
        ("forest", 4)
    ]

    # Create test data if it doesn't exist
    if !isdir(reference_path)
        println("Creating test dataset and reference results...")
        create_test_dataset()
        create_full_reference_dataset()
    else
        # Check if all reference files exist
        missing_refs = []
        for (tile, snfrac) in configurations
            ref_file = joinpath(reference_path, "reference_results_$(tile)_SNFRAC_$(snfrac).mat")
            if !isfile(ref_file)
                push!(missing_refs, (tile, snfrac))
            end
        end
        
        if !isempty(missing_refs)
            println("Missing reference files for configurations: $missing_refs")
            println("Creating missing reference results...")
            create_full_reference_dataset()
        end
    end
    
    # Run regression tests for each configuration
    for (tile, snfrac) in configurations
        @testset "$(tile) tile, SNFRAC=$(snfrac)" begin
            println("\n" * "="^60)
            println("Testing configuration: tile=$tile, SNFRAC=$snfrac")
            println("="^60)
            
            test_passed, max_differences = run_regression_test(test_data_path, tile, snfrac, tolerance=1e-10)
            
            @test test_passed            
        end
    end

end
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

    # Configuration matrix
    settings = [
        Dict(
            "tile" => "open",
            "config" => Dict("SNFRAC" => 0)
            ),
        Dict(
            "tile" => "open",
            "config" => Dict("SNFRAC" => 3)
            ),
        Dict(
            "tile" => "open",
            "config" => Dict("SNFRAC" => 4)
            ),
        Dict(
            "tile" => "forest",
            "config" => Dict("CANMOD" => 1, "EXCHNG" => 2, "SNFRAC" => 4, "ZOFFST" => 1),
            "params" => Dict("hfsn" => 0.3, "z0sn" => 0.01)
            )
    ]

    # Create test data if it doesn't exist
    if !isdir(reference_path)
        println("Creating test dataset and reference results...")
        create_test_dataset()
        create_full_reference_dataset()
    else
        # Check if all reference files exist
        missing_refs = []
        for setting in settings
            tile = setting["tile"]
            snfrac = setting["config"]["SNFRAC"]
            ref_file = joinpath(reference_path, "reference_results_$(tile)_SNFRAC_$(snfrac).mat")
            if !isfile(ref_file)
                push!(missing_refs, setting)
            end
        end
        
        if !isempty(missing_refs)
            println("Missing reference files for settings: $missing_refs")
            println("Creating missing reference results...")
            create_full_reference_dataset()
        end
    end
    
    # Run regression tests for each configuration
    for setting in settings
        tile = setting["tile"]
        snfrac = setting["config"]["SNFRAC"]
        @testset "$(tile) tile, SNFRAC=$(snfrac)" begin
            println("\n" * "="^60)
            println("Testing configuration: tile=$(setting["tile"]), snfrac=$(setting["config"]["SNFRAC"])")
            println("="^60)
            
            test_passed, max_differences = run_regression_test(test_data_path, setting, tolerance=1e-10)
            
            @test test_passed            
        end
    end

end
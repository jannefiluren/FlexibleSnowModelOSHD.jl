using Test
using FSMOSHD
using Serialization

# Absolute path to project folder
projdir = dirname(dirname(@__FILE__))

# Include regression test modules
include("regression_tests.jl")

# Test data paths
ref_file = joinpath(projdir, "test", "simulation_results.jls")

@testset "FSMOSHD Regression Tests" begin

    # Configuration matrix
    settings = [
        Dict(
            "tile" => "open",
            "config" => Dict("SNFRAC" => 0),
            "params" => Dict("wind_scaling" => 0.7)
            ),
        Dict(
            "tile" => "forest",
            "config" => Dict("CANMOD" => 1, "EXCHNG" => 2, "SNFRAC" => 4, "ZOFFST" => 1),
            "params" => Dict("hfsn" => 0.3, "z0sn" => 0.01, "wind_scaling" => 0.7)
            ),
        Dict(
            "tile" => "glacier",
            "config" => Dict("SNFRAC" => 0),
            "params" => Dict("wind_scaling" => 0.7)
        )
    ]

    # Create test data if it doesn't exist
    if !isfile(ref_file)
        println("Creating reference simulations and saving to file")
        simulation_refs = []
        for setting in settings
            push!(simulation_refs, run_simulations(setting))
        end
        serialize(ref_file, simulation_refs)
    else
        println("Loading reference simulations from file")
        simulation_refs = deserialize(ref_file)
    end

    # Run regression tests for each configuration
    tolerance = 1e-10
    for (setting, simulation_ref) in zip(settings, simulation_refs)

        test_passed = true

        tile = setting["tile"]

        println("Running regression test for tile = $tile")

        simulation_tst = run_simulations(setting)

        # Get list of variables to compare (exclude metadata)
        variables_to_compare = filter(k -> !(k in ["timestamps", "filenames"]), keys(simulation_tst))

        for variable in variables_to_compare
            
            # Check dimensions
            if size(simulation_tst[variable]) != size(simulation_ref[variable])
                println("ERROR: Dimension mismatch for $variable: current=$(size(simulation_tst[variable])), reference=$(size(simulation_ref[variable]))")
                test_passed = false
                continue
            end

            # Check values
            max_diff = maximum(abs.(simulation_tst[variable]-simulation_ref[variable]))
            if max_diff > tolerance
                println("FAIL for variable $variable: max_diff = $max_diff > $tolerance")
                test_passed = false
            end
        end
        
        if test_passed
            println("✅ Regression test passed for tile = $tile")
        else
            println("❌ Regression test failed for tile = $tile")
        end

        @test test_passed
        
    end

end
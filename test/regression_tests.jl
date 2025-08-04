using Test
using FSMOSHD
using MAT
using Dates


function load_reference_results(reference_path::String)
    reference_file = joinpath(reference_path, "reference_results.mat")
    
    if !isfile(reference_file)
        error("Reference results file not found: $reference_file")
    end
    
    reference_results = matread(reference_file)
    return reference_results
end


function run_regression_test(test_data_path::String; tolerance::Float64=1e-10)
    
    println("Running regression test...")
    println("Test data: $test_data_path")
    println("Tolerance: $tolerance")
    
    # Run current model to get simulation results
    println("Running current model...")
    current_results = run_snow_model(test_data_path)
    
    # Load reference results
    println("Loading reference results...")
    reference_path = joinpath(test_data_path, "reference")
    reference_results = load_reference_results(reference_path)
    
    # Compare results
    println("Comparing results...")
    test_passed = true
    max_differences = Dict{String, Float64}()
    
    # Get list of variables to compare (exclude metadata)
    variables_to_compare = filter(k -> !(k in ["timestamps", "filenames"]), keys(current_results))
    
    for variable in variables_to_compare
        if !haskey(reference_results, variable)
            println("WARNING: Variable $variable not found in reference results")
            continue
        end
        
        # Check that we have the same number of timesteps
        if length(current_results[variable]) != length(reference_results[variable])
            println("ERROR: Timestep mismatch for $variable: current=$(length(current_results[variable])), reference=$(length(reference_results[variable]))")
            test_passed = false
            continue
        end
        
        # Compare each timestep
        max_diff_var = 0.0
        for timestep in 1:length(current_results[variable])
            current_val = current_results[variable][timestep]
            reference_val = reference_results[variable][timestep]
            
            # Calculate maximum absolute difference for this timestep
            max_diff_timestep = maximum(abs.(current_val .- reference_val))
            max_diff_var = max(max_diff_var, max_diff_timestep)
            
            # Check if difference exceeds tolerance
            if max_diff_timestep > tolerance
                println("FAIL at timestep $timestep, variable $variable: max_diff = $max_diff_timestep > $tolerance")
                test_passed = false
            end
        end
        
        max_differences[variable] = max_diff_var
    end
    
    # Print summary
    println("\nRegression test summary:")
    println("=" ^ 50)
    for variable in sort(collect(keys(max_differences)))
        max_diff = max_differences[variable]
        status = max_diff <= tolerance ? "PASS" : "FAIL"
        println("$variable: max_diff = $(round(max_diff, sigdigits=3)) [$status]")
    end
    println("=" ^ 50)
    
    if test_passed
        println("✅ All regression tests PASSED")
    else
        println("❌ Some regression tests FAILED")
    end
    
    return test_passed, max_differences
end

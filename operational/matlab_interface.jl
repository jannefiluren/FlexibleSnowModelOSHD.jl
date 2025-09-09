#!/usr/bin/env julia

"""
MATLAB Interface for FSMOSHD Operational Runs

Command line interface to run_operational function for calling from MATLAB.

Usage:
    julia matlab_interface.jl <start_time> <end_time> [restart]

Arguments:
    start_time    Start time in ISO format (YYYY-MM-DDTHH:MM:SS)
    end_time      End time in ISO format (YYYY-MM-DDTHH:MM:SS)  
    restart       Optional: "true" to restart from existing state, "false" or omitted for new run

Examples:
    julia matlab_interface.jl 2024-09-01T06:00:00 2025-06-01T06:00:00
    julia matlab_interface.jl 2024-09-01T06:00:00 2025-06-01T06:00:00 true
"""

using FSMOSHD
using Dates
using ProgressMeter
using MAT

# Include the operational script
include("test_run_grid_operational.jl")

function main()
    # Check command line arguments
    if length(ARGS) < 2
        println("Error: Insufficient arguments")
        exit(1)
    end
    
    try
        # Parse start time
        start_str = ARGS[1]
        tstart = DateTime(start_str, "yyyy-mm-ddTHH:MM:SS")
        
        # Parse end time
        end_str = ARGS[2]
        tend = DateTime(end_str, "yyyy-mm-ddTHH:MM:SS")
        
        # Parse restart flag (optional, defaults to false)
        restart = false
        if length(ARGS) >= 3
            restart_str = lowercase(strip(ARGS[3]))
            if restart_str in ["true", "1", "yes"]
                restart = true
            elseif restart_str in ["false", "0", "no"]
                restart = false
            else
                restart = false  # Default to false for invalid input
            end
        end
        
        # Validate time range
        if tend <= tstart
            println("Error: End time must be after start time")
            exit(1)
        end
        
        # Call the operational function
        run_operational(tstart, tend, restart)
        
        # Simple success message
        println("Success: FSMOSHD run completed")
        exit(0)
        
    catch
        println("Error: FSMOSHD run failed")
        exit(1)
    end
end

# Run main function if script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
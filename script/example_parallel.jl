#!/usr/bin/env julia

"""
Example script showing how to use the parallel grid simulation.

Usage:
    julia --project=. -t 4 script/example_parallel.jl
    julia --project=. -t 8 script/example_parallel.jl
"""

using Pkg
Pkg.activate(".")

using Dates
using FSMOSHD
include("run_grid_simulation_parallel.jl")

println("="^50)
println("FSMOSHD Parallel Simulation Example")
println("="^50)
println("Available threads: $(Threads.nthreads())")
println()

# Example 1: Basic usage with default settings
println("Example 1: Basic parallel simulation")
run_grid_simulation_parallel(
    settings=Dict("tile" => "open"),
    times=DateTime(2024, 12, 1, 6):Hour(1):DateTime(2024, 12, 1, 9),
    verbose=true
)

# Example 2: Custom chunk count  
println("\nExample 2: Custom chunk count (half the threads)")
run_grid_simulation_parallel(
    settings=Dict("tile" => "open"), 
    times=DateTime(2024, 12, 1, 10):Hour(1):DateTime(2024, 12, 1, 12),
    nchunks=max(1, Threads.nthreads() ÷ 2),
    verbose=true
)

println("\n✅ Examples completed successfully!")
println("Check output in D:/julia/FSM_HS_julia/ directory")
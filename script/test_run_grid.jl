using Distributed

# Add worker processes (adjust based on available cores)
# Use nprocs() - 1 to leave one core for the main process
n_workers = min(4, Sys.CPU_THREADS - 1)  # Limit to 4 workers or available cores
if nprocs() == 1  # Only add workers if not already added
    addprocs(n_workers)
    println("Added $n_workers worker processes")
end

# Load required packages on all workers
@everywhere begin
    using ProgressMeter, MAT, Dates, FSMOSHD
    include("run_grid_simulation.jl")
end

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

println("Running $(length(settings)) configurations in parallel on $(nprocs()) processes")

# Run configurations in parallel
pmap(settings) do (setting)
    worker_id = myid()
    
    tile = setting["tile"]
    snfrac = setting["config"]["SNFRAC"]
    subfolder = uppercase(tile) * "_SNFRAC_" * string(snfrac)

    println("Worker $worker_id: Starting $tile SNFRAC=$snfrac -> $subfolder")

    start_time = time()
    run_grid_simulation(settings=setting, subfolder=subfolder, verbose=false)
    elapsed = time() - start_time

    println("Worker $worker_id: Completed $tile SNFRAC=$snfrac in $(round(elapsed/60, digits=1)) minutes")
end

println("All configurations completed")

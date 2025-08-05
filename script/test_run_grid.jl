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

configurations = [
    ("open", 0, "SNFRAC_0"),
    ("open", 3, "SNFRAC_3"),
    ("open", 4, "SNFRAC_4"),
    ("forest", 4, "FOREST"),
]

println("Running $(length(configurations)) configurations in parallel on $(nprocs()) processes")

# Run configurations in parallel
pmap(configurations) do (tile, snfrac, subfolder)
    worker_id = myid()
    println("Worker $worker_id: Starting $tile SNFRAC=$snfrac -> $subfolder")
    
    start_time = time()
    run_grid_simulation(tile=tile, snfrac=snfrac, subfolder=subfolder, verbose=false)
    elapsed = time() - start_time
    
    println("Worker $worker_id: Completed $tile SNFRAC=$snfrac in $(round(elapsed/60, digits=1)) minutes")
end

println("All configurations completed")

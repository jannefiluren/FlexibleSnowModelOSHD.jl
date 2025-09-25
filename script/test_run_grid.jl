using Distributed

# Configure simulations
settings_common = Dict(
    "met_type" => :ICON,
    "met_folder" => "K:/DATA_ICON/OUTPUT_OSHD_0250/PROCESSED_ANALYSIS/ICON_1EFA",
    "met_prefix" => "ICON",
    "tvt_folder" => "K:/OSHD_AUX/DATA_CANRAD/OUTPUT_OSHD_0250/CR_2410_static",
    "lus_file" => "K:/OSHD_AUX/DATA_LUS/OSHD_LUS_0250.mat",
    "output_vars" => ["snowdepth", "fsnow", "SWE", "Sliq_out", "Tsrf", "Melt", "Roff", "meltflux_out"],
    )

settings_specific = [
    Dict(
        "tile" => "open",
        "config" => Dict("SNFRAC" => 0),
        "params" => Dict("wind_scaling" => 0.7),
        "out_folder" => "D:/MODEL_DATA_FSM/FSM_HS_julia/OPEN_SNFRAC_3",
    ),
    Dict(
        "tile" => "open",
        "config" => Dict("SNFRAC" => 0),
        "params" => Dict("wind_scaling" => 0.7, "scf_a" => 0.3193*250^0.103, "scf_b"=>0.5333*250^0.0389),
        "out_folder" => "D:/MODEL_DATA_FSM/FSM_HS_julia/OPEN_SNFRAC_4",
    ),
    Dict(
        "tile" => "open",
        "config" => Dict("SNFRAC" => 0),
        "params" => Dict("wind_scaling" => 0.7, "scf_a" => 0.3193*100^0.103, "scf_b"=>0.5333*100^0.0389),
        "out_folder" => "D:/MODEL_DATA_FSM/FSM_HS_julia/OPEN_SNFRAC_5",
    ),
]

settings = [merge(settings_common, settings) for settings in settings_specific]

# Add worker processes (adjust based on available cores)
# Use nprocs() - 1 to leave one core for the main process
n_workers = min(length(settings), Sys.CPU_THREADS - 1)  # Limit to length(settings) workers or available cores
if nprocs() == 1  # Only add workers if not already added
    addprocs(n_workers)
    println("Added $n_workers worker processes")
end

# Load required packages on all workers
@everywhere begin
    using ProgressMeter, MAT, Dates, FSMOSHD
    include("run_grid_simulation.jl")
end

# Run configurations in parallel
println("Running $(length(settings)) configurations in parallel on $(nprocs()) processes")

pmap(settings) do (setting)
    worker_id = myid()
    
    tile = setting["tile"]
    snfrac = setting["config"]["SNFRAC"]
    subfolder = splitpath(setting["out_folder"])[end]

    println("Worker $worker_id: Starting $tile SNFRAC=$snfrac -> $subfolder")

    start_time = time()
    run_grid_simulation(settings=setting, verbose=false)
    elapsed = time() - start_time

    println("Worker $worker_id: Completed $tile SNFRAC=$snfrac in $(round(elapsed/60, digits=1)) minutes")
end

println("All configurations completed")

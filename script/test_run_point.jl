using Dates

include("run_point_simulation.jl")

# Run the point simulation with default settings
results = run_point_simulation(
    times = DateTime(2024, 9, 1, 6):Hour(1):DateTime(2025, 06, 12, 6),
    landuse_file = "K:/OSHD_AUX/DATA_LUS/OSHD_LUS_STAT.mat",
    meteo_base_path = "K:/DATA_ICON/OUTPUT_OSHD_STAT/PROCESSED_ANALYSIS/ICON_1EFA",
    reference_path = "D:/julia/FSM_HS_matlab/LATEST_00h_RUN/OUTPUT_OSHD_STAT/RESULTS_01h_opn",
    base_folder = "D:/julia",
    write_outputs = true,
    verbose = true
)

println("Simulation completed successfully!")
println("Snow depth dimensions: ", size(results.snowdepth))
println("Reference dimensions: ", size(results.snowdepth_ref))

# Compute maximum difference between simulations
max_diff = compute_max_difference(results.snowdepth, results.snowdepth_ref)
println("Maximum difference: ", round(max_diff, digits=8), " m")
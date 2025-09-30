using FSMOSHD, Dates


println("------------- OSHD DEFAULT RUN -------------")

times = DateTime(2024, 12, 1, 6):Hour(1):DateTime(2024, 12, 2, 6)

run_grid_simulation(times=times)


println("------------- GLAMOS RUN -------------")

settings = Dict(
    "tile" => "open",
    "config" => Dict("SNFRAC" => 0),
    "params" => Dict("wind_scaling" => 1.0),
    "met_type" => :GLAMOS,
    "met_folder" => "D:/julia/glamos/metdata",
    "met_prefix" => "COSMO",
    "tvt_folder" => "K:/OSHD_AUX/DATA_CANRAD/OUTPUT_OSHD_0250/CR_2410_static",
    "lus_file" => "D:/julia/glamos/lusdata/GPMO_LUS_2025_25.mat",
    "out_folder" => string(@__DIR__, "/../../FSM_HS_julia/default"),
    "output_vars" => ["snowdepth", "fsnow", "Roff", "meltflux_out"],
    )

times = DateTime(2022, 12, 1, 6):Hour(1):DateTime(2022, 12, 2, 6)

run_grid_simulation(settings=settings, times=times)


println("------------- SNOWINFLOW RUN -------------")

settings = Dict(
    "tile" => "open",
    "config" => Dict("SNFRAC" => 0),
    "params" => Dict("wind_scaling" => 1.0),
    "met_type" => :COSMO,
    "met_folder" => "D:/snowinflow_project/snowinflow_data/model_data/fsm_data/DATA_MET_NORDIC/OUTPUT_ULLA_FORRE_0900/PROCESSED_ANALYSIS/",
    "met_prefix" => "METNORDIC",
    "lus_file" => "D:/snowinflow_project/snowinflow_data/model_data/fsm_data/SNOWINFLOW_AUX/DATA_LUS/ULLA_FORRE_LUS_0900.mat",
    "out_folder" => string(@__DIR__, "/../../FSM_HS_julia/default"),
    "output_vars" => ["snowdepth", "fsnow", "Roff", "meltflux_out"],
    )

times = DateTime(2022, 12, 1, 6):Hour(1):DateTime(2022, 12, 2, 6)

run_grid_simulation(settings=settings, times=times)

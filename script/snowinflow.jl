include("run_grid_simulation.jl")

lus_file = "D:/snowinflow_project/snowinflow_data/model_data/fsm_data/SNOWINFLOW_AUX/DATA_LUS/ULLA_FORRE_LUS_0900.mat"
met_folder = "D:/snowinflow_project/snowinflow_data/model_data/fsm_data/DATA_MET_NORDIC/OUTPUT_ULLA_FORRE_0900"
subfolder = "model_data/fsm_data/FSM_julia"

times = DateTime(2021, 9, 1, 6):Hour(1):DateTime(2022, 9, 1, 6)

settings = Dict(
  "tile" => "open",
  "config" => Dict(
    "OSHDTN" => 0,
    "SNFRAC" => 0,
  ),
  "params" => Dict(
    "zT" => 2.0,
  )
)

run_grid_simulation(; lus_file=lus_file, met_folder=met_folder, times=times, settings=settings, subfolder=subfolder)
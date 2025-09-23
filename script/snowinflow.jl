include("run_grid_simulation.jl")
include("user_paths.jl")

using GeoIO
using Dates

resolution = "0900"
region = "ULLA_FORRE"


lus_file = local_folder * "/model_data/fsm_data/SNOWINFLOW_AUX/DATA_LUS/$(region)_LUS_$(resolution).mat"
met_folder = local_folder * "/model_data/fsm_data/DATA_MET_NORDIC/OUTPUT_$(region)_$(resolution)"
out_folder = local_folder * "/model_data/fsm_data/FSM_julia/_$(resolution)"

lidar_file = "C:/Users/opreljen/OneDrive - Eidg. Forschungsanstalt WSL/Statkraft/lidar_flights/2022/data/masked/resampled/ASO_Kvilldal_Mosaic_2022Mar19_snowdepth_3m_maskedLakes_resampled_900.tif"
lidar_import = GeoIO.load(lidar_file) #for some reason CRS is not in this struct.
lidar_data = reshape(lidar_import.color, lidar_import.geometry.topology.dims)'
lidar_time = DateTime(2022, 3, 19, 12)

times = DateTime(2021, 9, 1, 6):Hour(1):DateTime(2022, 4, 1, 6)

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


scaling_grid = ones(size(lidar_data))     


#put scaling grid in settings
scaling_settings = Dict(
    "sfscaling" => scaling_grid, #TODO: Set default somewhere
)
merge!(settings, scaling_settings)

fsm_hs = run_grid_simulation(; lus_file=lus_file, met_folder=met_folder, times=times, 
                            settings=settings, out_folder=out_folder, hs_out_times=[lidar_time])
fsm_hs = dropdims(fsm_hs, dims=3) #remove singleton dimension    

scaling_grid .= ifelse.(isnan.(lidar_data), ones(size(fsm_hs)), lidar_data ./ fsm_hs) .* scaling_grid #update scaling grid
scaling_settings = Dict(
    "sfscaling" => scaling_grid, #TODO: Set default somewhere
)
merge!(settings, scaling_settings)


using Plots
plotlyjs()  # or gr(), depending on your backend preference

display(heatmap(scaling_grid; color=:viridis, aspect_ratio=:equal, title="Scaling Grid", colorbar=true))
display(heatmap(fsm_hs; color=:viridis, aspect_ratio=:equal, title="FSM hs", colorbar=true))

fsm_hs = run_grid_simulation(; lus_file=lus_file, met_folder=met_folder, times=times, 
                            settings=settings, out_folder=out_folder, hs_out_times=[lidar_time])
fsm_hs = dropdims(fsm_hs, dims=3) #remove singleton dimension   

display(heatmap(fsm_hs; color=:viridis, aspect_ratio=:equal, title="FSM hs", colorbar=true))

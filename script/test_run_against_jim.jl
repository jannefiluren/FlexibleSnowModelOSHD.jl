using FSMOSHD
using Dates
using MAT
using Printf

# Helper functions

include("test_run_against_jim_functions.jl")

# Settings

config = Dict()
config["base_folder"] = "D:/julia"
config["test_setup"] = false
config["test_drive"] = false
config["test_radiation"] = false
config["test_thermal"] = false
config["test_sfexch"] = false
config["test_ebalsrf"] = false
config["test_snow"] = false
config["test_soil"] = false
config["test_snow_final"] = true

time_initialize = DateTime(2024,9,1,6,0,0)
time_start = DateTime(2024,9,1,6,0,0)
time_end = DateTime(2025,6,1,6,0,0)

station = "MCH.JUN2"
matfiles = false

# Compile model

compile_fortran_code(config["base_folder"])

# Prepare landuse

landuse = prepare_landuse_stations(station)

# Loop over time

times = time_start:Dates.Day(1):time_end

for currtime in times

  global fsm, meteo

  # Run model from matlab
  
  if currtime == time_initialize
    cmd = `matlab -batch "cd('D:\julia\FSMOSHD\script'); test_run_point_function($(year(currtime)), $(month(currtime)), $(day(currtime)), 'initialize', '$(station)')"`
  else
    cmd = `matlab -batch "cd('D:\julia\FSMOSHD\script'); test_run_point_function($(year(currtime)), $(month(currtime)), $(day(currtime)), 'reinitialize', '$(station)')"`
  end
  
  pr = run(cmd; wait=true)

  # Setup

  if currtime == times[1]
    if matfiles
      fsm = setup_matfiles(Float32, Int32, landuse, sum(landuse["is_domain"]), 1)
    else
      fsm = setup_binfiles(Float32, Int32, joinpath(config["base_folder"], "FSM_HS_single/bin_files"))
    end
    meteo = MET{Float32,Int32}(Nx=fsm.Nx, Ny=fsm.Ny)
  end
  
  # Run tests
  
  failure, meteo = test_run_against_jim_binfiles(config, fsm, meteo, currtime, station, matfiles)

  if failure
    @error "Mismatch between results"
    break
  end

end






# using FSMOSHD
# using Dates

# # Helper functions

# include("test_run_against_jim_functions.jl")

# # Settings

# config = Dict()
# config["base_folder"] = "D:/julia"
# config["test_setup"] = false
# config["test_drive"] = false
# config["test_radiation"] = false
# config["test_thermal"] = false
# config["test_sfexch"] = true
# config["test_ebalsrf"] = false
# config["test_snow"] = false
# config["test_soil"] = false
# config["test_snow_final"] = false

# # Run tests

# test_run_against_jim(config)
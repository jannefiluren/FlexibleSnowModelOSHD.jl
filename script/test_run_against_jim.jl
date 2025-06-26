using FSMOSHD
using Dates

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

# Compile model

compile_fortran_code(config["base_folder"])

# Loop over time

time_initialize = Date(2024,9,1)
time_start = Date(2024,9,1)
time_end = Date(2025,6,1)

times = time_start:time_end

for time in times

  # Run model from matlab
  
  if time == time_initialize
    cmd = `matlab -batch "cd('D:\julia\FSMOSHD\script'); test_run_point_function($(year(time)), $(month(time)), $(day(time)), 'initialize')"`
  else
    cmd = `matlab -batch "cd('D:\julia\FSMOSHD\script'); test_run_point_function($(year(time)), $(month(time)), $(day(time)), 'reinitialize')"`
  end
  
  pr = run(cmd; wait=true)
  
  # Run tests
  
  failure = test_run_against_jim_binfiles(config)

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
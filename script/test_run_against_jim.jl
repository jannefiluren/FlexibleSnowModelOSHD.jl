using FSMOSHD
using Dates

# Helper functions

include("test_functions.jl")


# Settings

base_folder = "D:/julia"
test_setup = false
test_drive = false
test_radiation = false
test_thermal = false
test_sfexch = false
test_ebalsrf = false
test_snow = true
test_soil = false


# Compile and run original fortran code

include("test_run_against_jim_compile_run.jl")


# Create object from bin files

fsm = setup_original(Float32, Int32, joinpath(base_folder, "FSM_HS/bin_files"));
meteo = MET{Float32, Int32}(Nx=fsm.Nx, Ny=fsm.Ny)

# Open files

io = open_files(base_folder, fsm)

# Test setup

if test_setup
  verify_setup(base_folder, fsm)
end

# Loop over time

for i = 1:24

  # Read drive

  drive_original!(io, meteo, fsm)

  if test_drive
    verify_drive(base_folder, meteo)
  end

  # Run radiation

  t = DateTime(meteo.year[1], meteo.month[1], meteo.day[1], meteo.hour[1], 00, 00)

  radiation(fsm, meteo, t)

  if test_radiation
    verify_radiation(base_folder, fsm)
  end

  # Run thermal

  thermal(fsm)

  if test_thermal
    verify_thermal(base_folder, fsm)
  end

  # Run sfexch and ebalsrf

  for i in 1:fsm.Nitr
    sfexch(fsm, meteo)
    ebalsrf(fsm, meteo)
  end

  if test_sfexch
    verify_sfexch(base_folder, fsm)
  end

  if test_ebalsrf
    verify_ebalsrf(base_folder, fsm)
  end

  # Run snow

  snow(fsm, meteo, t)

  if test_snow
    verify_snow(base_folder, fsm)
  end

  # Run soil

  soil(fsm)

  if test_soil
    verify_soil(base_folder, fsm)
  end

end
using ProgressMeter
using MAT
using Dates
using FSMOSHD
using CSV
using Tables

# Settings

SNFRAC = 0

base_folder = "D:/julia"
subfolder = "SNFRAC_" * string(SNFRAC)

# Helper functions

searchdir(path, key) = filter(x -> occursin(key, x), readdir(path))

# Read landuse data

landuse = prepare_landuse_grid()


# Setup model

Nx = size(landuse["dem"]["data"], 1)
Ny = size(landuse["dem"]["data"], 2)

fsm = setup_matfiles(Float32, Int32, landuse, Nx, Ny, SNFRAC = SNFRAC)
met_curr = MET{Float32, Int32}(Nx=Nx, Ny=Ny)

# Run model

mkpath(joinpath(base_folder, "FSM_HS_julia", subfolder))

times = DateTime(2024, 9, 1, 6):Hour(1):DateTime(2025, 06, 12, 6)

for (i, t) in enumerate(times)
# @showprogress "Running snow model..." for (i, t) in enumerate(times)

  @time begin
  
    drive_grid!(met_curr, fsm, t)
  
    radiation(fsm, met_curr, t)
  
    thermal(fsm)
  
    for i in 1:fsm.Nitr
      sfexch(fsm, met_curr)
      ebalsrf(fsm, met_curr)
    end
  
    snow(fsm, met_curr, t)
  
    soil(fsm)

    if hour(t) == 5
    
      matwrite(joinpath(base_folder, "FSM_HS_julia", subfolder, Dates.format(t + Dates.Hour(1), "yyyymmddHHMM") * "_output.mat"),
      Dict(
        "hs" => dropdims(sum(fsm.Ds, dims=1), dims=1).*fsm.fsnow,
        "fsnow" => fsm.fsnow
        ); compress = true)
        
    end
    
  end
  
end

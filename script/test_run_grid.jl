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

fsm = setup_matfiles_grid(Float32, Int32, landuse, Nx, Ny, SNFRAC = SNFRAC)
met_curr = MET{Float32, Int32}(Nx=Nx, Ny=Ny)

Sf24h = zeros(size(met_curr.Sf24h))
Sf_history = zeros(size(met_curr.Sf_history))

# Run model

mkpath(joinpath(base_folder, "FSM_HS_julia", subfolder))

times = DateTime(2024, 9, 1, 6):Hour(1):DateTime(2025, 06, 12, 6)

for (i, t) in enumerate(times)

  @time begin

    # Prepare forcing data

    folder = joinpath("K:/DATA_ICON/OUTPUT_OSHD_0250/PROCESSED_ANALYSIS/ICON_1EFA", Dates.format(t, "yyyy.mm"))
    filename = searchdir(folder, "ICONDATA_" * Dates.format(t, "yyyymmddHHMM") * "_C1EFA_")

    met_single = matread(joinpath(folder, filename[1]))

    Sdir = met_single["sdri"]["data"]     # "direct shortwave radiation, per inclined surface area, within topography, above canopy"
    Sdif = met_single["sdfd"]["data"]     # "diffuse shortwave radiation, per horizontal surface area, within topography, above canopy"
    Sdird = met_single["sdrd"]["data"]    # "direct shortwave radiation, per horizontal surface area, within topography, above canopy"
    LW = met_single["lwrs"]["data"]       # "longwave radiation, above topography"
    Rf = met_single["prfc"]["data"]       # "rainfall"
    Sf = met_single["psfc"]["data"]       # "snowfall (snow + graupel)"
    Ta = met_single["tais"]["data"]       # "air temperature"
    RH = met_single["rhus"]["data"]       # "relative humidity"
    Ua = met_single["wnss"]["data"]       # "wind speed"
    Ps = met_single["pail"]["data"]       # "local air pressure"

    met_curr.Sdir[:, :] .= Sdir
    met_curr.Sdif[:, :] .= Sdif
    met_curr.Sdird[:, :] = Sdird
    met_curr.LW[:, :] .= LW
    met_curr.Sf[:, :] .= Sf
    met_curr.Rf[:, :] .= Rf
    met_curr.Ta[:, :] .= Ta
    met_curr.RH[:, :] .= RH
    met_curr.Ua[:, :] .= Ua
    met_curr.Ps[:, :] .= Ps

    curr_hour = Dates.value(Hour(t)) + 1
    Sf24h .+= Sf
    Sf24h .-= Sf_history[:,:,curr_hour]
    Sf_history[:,:,curr_hour] = Sf
    met_curr.Sf24h[:, :] .= Sf24h

    # Run model
    
    drive!(fsm, met_curr)
  
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

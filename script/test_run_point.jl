using ProgressMeter
using MAT
using Dates
using FSMOSHD
using CSV
using Tables


# Settings

base_folder = "D:/julia"


# Helper functions

searchdir(path, key) = filter(x -> occursin(key, x), readdir(path))


# Read landuse data

landuse = matread("K:/OSHD_AUX/DATA_LUS/OSHD_LUS_STAT.mat")
nstat = length(landuse["acro"])


# Read input data

times = DateTime(2024, 9, 1, 6):Hour(1):DateTime(2025, 06, 12, 6)

Sdir = zeros(length(times), nstat)
Sdif = zeros(length(times), nstat)
Sdird = zeros(length(times), nstat)
LW = zeros(length(times), nstat)
Sf = zeros(length(times), nstat)
Rf = zeros(length(times), nstat)
Ta = zeros(length(times), nstat)
RH = zeros(length(times), nstat)
Ua = zeros(length(times), nstat)
Ps = zeros(length(times), nstat)

@showprogress "Reading meteo data..." for (i, t) in enumerate(times)
    folder = joinpath("K:/DATA_ICON/OUTPUT_OSHD_STAT/PROCESSED_ANALYSIS/ICON_1EFA", Dates.format(t, "yyyy.mm"))
    filename = searchdir(folder, "ICONDATA_" * Dates.format(t, "yyyymmddHHMM") * "_C1EFA_")

    met_single = matread(joinpath(folder, filename[1]))

    Sdir[i, :] = met_single["sdri"]["data"]     # "direct shortwave radiation, per inclined surface area, within topography, above canopy"
    Sdif[i, :] = met_single["sdfd"]["data"]     # "diffuse shortwave radiation, per horizontal surface area, within topography, above canopy"
    Sdird[i, :] = met_single["sdrd"]["data"]    # "direct shortwave radiation, per horizontal surface area, within topography, above canopy"
    LW[i, :] = met_single["lwrs"]["data"]       # "longwave radiation, above topography"
    Rf[i, :] = met_single["prfc"]["data"]       # "rainfall"
    Sf[i, :] = met_single["psfc"]["data"]       # "snowfall (snow + graupel)"
    Ta[i, :] = met_single["tais"]["data"]       # "air temperature"
    RH[i, :] = met_single["rhus"]["data"]       # "relative humidity"
    Ua[i, :] = met_single["wnss"]["data"]       # "wind speed"
    Ps[i, :] = met_single["pail"]["data"]       # "local air pressure"
end


# Prepare landuse

landuse = prepare_landuse("K:/OSHD_AUX/DATA_LUS/OSHD_LUS_STAT.mat")
landuse_cropped = crop_landuse_to_domain(landuse)

# Setup model
nstat_domain = sum(landuse["is_domain"]["data"])
fsm = setup(Float32, Int32, landuse_cropped, nstat_domain, 1)
met_curr = MET{Float32, Int32}(Nx=nstat)

Sf24h = zeros(size(met_curr.Sf24h))
Sf_history = zeros(size(met_curr.Sf_history))

# Run model

snowdepth = zeros(length(times), nstat)

@showprogress "Running snow model..." for (i, t) in enumerate(times)

  # Forcing data

  met_curr.Sdir[:, :] .= Sdir[i, :]
  met_curr.Sdif[:, :] .= Sdif[i, :]
  met_curr.Sdird[:, :] = Sdird[i, :]
  met_curr.LW[:, :] .= LW[i, :]
  met_curr.Sf[:, :] .= Sf[i, :]
  met_curr.Rf[:, :] .= Rf[i, :]
  met_curr.Ta[:, :] .= Ta[i, :]
  met_curr.RH[:, :] .= RH[i, :]
  met_curr.Ua[:, :] .= Ua[i, :]
  met_curr.Ps[:, :] .= Ps[i, :]

  curr_hour = Dates.value(Hour(t)) + 1
  Sf24h .+= Sf[i, :]
  Sf24h .-= Sf_history[:,:,curr_hour]
  Sf_history[:,:,curr_hour] = Sf[i, :]
  met_curr.Sf24h[:, :] .= Sf24h
  
  # Run model

  drive!(fsm, met_curr)

  radiation!(fsm, met_curr, t)

  thermal!(fsm)

  for i in 1:fsm.Nitr
    sfexch!(fsm, met_curr)
    ebalsrf!(fsm, met_curr)
  end

  snow!(fsm, met_curr, t)

  soil!(fsm)

  # Output data

  snowdepth[i,:] = dropdims(sum(fsm.Ds, dims=1), dims=3)

end


# Load reference results

snowdepth_ref = similar(snowdepth)

@showprogress "Reading reference results..." for (i, t) in enumerate(times)

  folder = "D:/julia/FSM_HS_all/LATEST_00h_RUN/OUTPUT_OSHD_STAT/RESULTS_01h_opn"
  filename = searchdir(folder, "MODELDATA_" * Dates.format(t, "yyyymmddHHMM"))

  file = matopen(joinpath(folder, filename[1]))
  hsnt = read(file, "hsnt") 
  close(file)

  snowdepth_ref[i, :] = hsnt["data"]

end


# Write results to files

CSV.write(joinpath(base_folder,"FSM_HS_all","snowdepth_julia.csv"), Tables.table(snowdepth))
CSV.write(joinpath(base_folder,"FSM_HS_all","snowdepth_matlab.csv"), Tables.table(snowdepth_ref))
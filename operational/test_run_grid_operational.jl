using ProgressMeter
using MAT
using Dates
using FSMOSHD


# tstart = DateTime(2024,9,1,6)
# tend = DateTime(2025,6,1,6)
# run_operational(tstart, tend)


function run_operational(tstart::DateTime, tend::DateTime, restart::Bool=false)

  # Settings

  base_folder = "D:/julia/FSM_HS_julia_operational"

  settings = Dict("tile" => "open", "config" => Dict("SNFRAC" => 0), "params" => Dict("wind_scaling" => 0.7))
  
  # Read landuse data

  landuse = prepare_landuse("K:/OSHD_AUX/DATA_LUS/OSHD_LUS_0250.mat")

  # Setup model

  Nx, Ny = size(landuse["dem"]["data"])

  fsm = setup(Float32, Int32, landuse, Nx, Ny, settings)
  met = MET{Float32, Int32}(Nx=Nx, Ny=Ny)

  Sf24h = zeros(size(met.Sf24h))
  Sf_history = zeros(size(met.Sf_history))

  if restart
    read_states!(fsm, Sf24h, Sf_history, base_folder, tstart)
  end

  # Create states and results directory

  mkpath(joinpath(base_folder, "results"))
  mkpath(joinpath(base_folder, "states"))

  # Run model

  times = tstart:Hour(1):tend

  @showprogress "Running snow model..." for t in times

    # Prepare forcing data

    folder_icon = joinpath("K:/DATA_ICON/OUTPUT_OSHD_0250/PROCESSED_ANALYSIS/ICON_1EFA", Dates.format(t, "yyyy.mm"))
    filename_icon = searchdir(folder_icon, "ICONDATA_" * Dates.format(t, "yyyymmddHHMM"))
    met_single = matread(joinpath(folder_icon, filename_icon[1]))

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

    met.Sdir[:, :] .= Sdir
    met.Sdif[:, :] .= Sdif
    met.Sdird[:, :] = Sdird
    met.LW[:, :] .= LW
    met.Sf[:, :] .= Sf
    met.Rf[:, :] .= Rf
    met.Ta[:, :] .= Ta
    met.RH[:, :] .= RH
    met.Ua[:, :] .= Ua
    met.Ps[:, :] .= Ps

    curr_hour = Dates.value(Hour(t)) + 1
    Sf24h .+= Sf
    Sf24h .-= Sf_history[:, :, curr_hour]
    Sf_history[:, :, curr_hour] = Sf
    met.Sf24h[:, :] .= Sf24h

    # Forest-specific TVT data loading

    if settings["tile"] == "forest"
        folder_tvt = "K:/OSHD_AUX/DATA_CANRAD/OUTPUT_OSHD_0250/CR_2410_static/YYYY." * Dates.format(t, "mm")
        filename_tvt = searchdir(folder_tvt, "CANRAD_" * Dates.format(t, "mmddHHMM"))
        tvt_single = matread(joinpath(folder_tvt, filename_tvt[1]))
        met.Tv[:, :] .= tvt_single["stdx"]["data"]
    end

    # Run model for this timestep

    step!(fsm, met, t)

    # Write states and results

    if hour(t) == 5
      write_results(fsm, base_folder, t)
      write_states(fsm, Sf24h, Sf_history, base_folder, t)
    end

  end

end

function write_results(fsm, folder, t)

  matwrite(joinpath(folder, "results", Dates.format(t + Dates.Hour(1), "yyyymmddHHMM") * "_output.mat"),
        Dict(
          "hs" => dropdims(sum(fsm.Ds, dims=1), dims=1) .* fsm.fsnow,
          "fsnow" => fsm.fsnow
        ); compress=true)

end

function write_states(fsm, Sf24h, Sf_history, folder, t)

  variables = ["albs", "Ds", "Nsnow", "Qcan", "rgrn", "histowet", "Sice", "Sliq",
    "Sveg", "Tcan", "theta", "Tsnow", "Tsoil", "Tsrf", "fsnow", "Tveg", "snowdepthmin",
    "snowdepthmax", "snowdepthhist", "swemin", "swemax", "swehist", "fsky_terr"]

  states = Dict()
  for variable in variables
    state = Dict()
    state["data"] = getfield(fsm, Symbol(variable))
    state["desc"] = "state variable: " * variable
    states[variable] = state
  end

  state = Dict()
  state["data"] = Sf24h
  state["desc"] = "state variable: Sf24h"
  states["Sf24h"] = state

  state = Dict()
  state["data"] = Sf_history
  state["desc"] = "state variable: Sf_history"
  states["Sf_history"] = state

  file = joinpath(folder, "states", Dates.format(t + Dates.Hour(1), "yyyymmddHHMM") * "_states.mat")

  matwrite(file, states, ; compress=true)

end

function read_states!(fsm, Sf24h, Sf_history, folder, t)

  file = joinpath(folder, "states", Dates.format(t, "yyyymmddHHMM") * "_states.mat")

  states = matread(file)

  for variable in keys(states)
    if variable == "Sf24h"
      Sf24h .= states["Sf24h"]["data"]
    elseif variable == "Sf_history"
      Sf_history .= states["Sf_history"]["data"]
    else
      state = states[variable]["data"]
      setfield!(fsm, Symbol(variable), state)
    end
  end

end

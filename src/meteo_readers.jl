function read_meteo!(met::MET{Tf,Ti}, fsm::FSM{Tf,Ti}, t::DateTime, settings::Dict) where {Tf<:Real,Ti<:Integer}
  read_meteo!(met, fsm, t, settings, Val(settings["met_type"]))
end


function read_meteo!(met::MET{Tf,Ti}, fsm::FSM{Tf,Ti}, t::DateTime, settings::Dict, ::Val{:ICON}) where {Tf<:Real,Ti<:Integer}

  # Prepare forcing data
  folder_met = joinpath(settings["met_folder"], Dates.format(t, "yyyy.mm"))
  filename_met = searchdir(folder_met, settings["met_prefix"] * "DATA_" * Dates.format(t, "yyyymmddHHMM"))
  met_single = matread(joinpath(folder_met, filename_met[1]))

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

  if fsm.SNTRAN == 1
    if haskey(met_single, "wnsd")
      Ua = met_single["wnsd"]["data"]       # "wind speed"
    else
      Ua = met_single["wnss"]["data"]       # "wind speed"
    end
    if haskey(met_single, "wndd")
      Udir = met_single["wndd"]["data"]     # "wind direction"
    else
      Udir = met_single["wnds"]["data"]     # "wind direction"
    end
    met.Ua[:, :] .= Ua
    met.Udir[:, :] .= Udir
    Udir[isnan.(Udir)] .= 0
  end

  # Update 24-hour snowfall tracking
  curr_hour = Dates.value(Hour(t)) + 1
  met.Sf24h_f64 .+= Sf
  met.Sf24h_f64 .-= met.Sf_history_f64[:, :, curr_hour]
  met.Sf_history_f64[:, :, curr_hour] = Sf
  met.Sf24h[:, :] .= met.Sf24h_f64

  # Forest-specific time-varying transmissivity data loading
  if fsm.TILE == "forest"
    folder_tvt = joinpath(settings["tvt_folder"], "YYYY." * Dates.format(t, "mm"))
    filename_tvt = searchdir(folder_tvt, "CANRAD_" * Dates.format(t, "mmddHHMM"))
    tvt_single = matread(joinpath(folder_tvt, filename_tvt[1]))
    met.Tv[:, :] .= tvt_single["stdx"]["data"]
  end

end

function read_meteo!(met::MET{Tf,Ti}, fsm::FSM{Tf,Ti}, t::DateTime, settings::Dict, ::Val{:COSMO}) where {Tf<:Real,Ti<:Integer}

  # Prepare forcing data
  folder_met = joinpath(settings["met_folder"], Dates.format(t, "yyyy.mm"))
  filename_met = searchdir(folder_met, settings["met_prefix"] * "DATA_" * Dates.format(t, "yyyymmddHHMM"))
  met_single = matread(joinpath(folder_met, filename_met[1]))

  Sdir = met_single["sdrd"]["data"]     # "direct shortwave radiation, per inclined surface area, within topography, above canopy"
  Sdif = met_single["sdfd"]["data"]     # "diffuse shortwave radiation, per horizontal surface area, within topography, above canopy"
  LW = met_single["lwrs"]["data"]       # "longwave radiation, above topography"
  Ptot = met_single["prcs"]["data"]     # "precipitation"
  Ta = met_single["tais"]["data"]       # "air temperature"
  RH = met_single["rhus"]["data"]       # "relative humidity"
  Ua = met_single["wnss"]["data"]       # "wind speed"
  Ps = met_single["pail"]["data"]       # "local air pressure"

  thres_prec = 1.04 + 273.15
  m_prec = 0.15
  Sf = 1.0 ./ (1 .+ exp.((Ta .- thres_prec) ./ m_prec)) .* Ptot
  Rf = Ptot .- Sf

  met.Sdir[:, :] .= Sdir
  met.Sdif[:, :] .= Sdif
  met.LW[:, :] .= LW
  met.Sf[:, :] .= Sf
  met.Rf[:, :] .= Rf
  met.Ta[:, :] .= Ta
  met.RH[:, :] .= RH
  met.Ua[:, :] .= Ua
  met.Ps[:, :] .= Ps

  if fsm.SNTRAN == 1
    if haskey(met_single, "wnsd")
      Ua = met_single["wnsd"]["data"]       # "wind speed"
    else
      Ua = met_single["wnss"]["data"]       # "wind speed"
    end
    if haskey(met_single, "wndd")
      Udir = met_single["wndd"]["data"]     # "wind direction"
    else
      Udir = met_single["wnds"]["data"]     # "wind direction"
    end
    met.Ua[:, :] .= Ua
    met.Udir[:, :] .= Udir
    Udir[isnan.(Udir)] .= 0
  end

  # Update 24-hour snowfall tracking
  curr_hour = Dates.value(Hour(t)) + 1
  met.Sf24h_f64 .+= Sf
  met.Sf24h_f64 .-= met.Sf_history_f64[:, :, curr_hour]
  met.Sf_history_f64[:, :, curr_hour] = Sf
  met.Sf24h[:, :] .= met.Sf24h_f64

  # Forest-specific time-varying transmissivity data loading
  if fsm.TILE == "forest"
    folder_tvt = joinpath(settings["tvt_folder"], "YYYY." * Dates.format(t, "mm"))
    filename_tvt = searchdir(folder_tvt, "CANRAD_" * Dates.format(t, "mmddHHMM"))
    tvt_single = matread(joinpath(folder_tvt, filename_tvt[1]))
    met.Tv[:, :] .= tvt_single["stdx"]["data"]
  end

end

function read_meteo!(met::MET{Tf,Ti}, fsm::FSM{Tf,Ti}, t::DateTime, settings::Dict, ::Val{:GLAMOS}) where {Tf<:Real,Ti<:Integer}

  # Prepare forcing data
  folder_met = joinpath(settings["met_folder"], Dates.format(t, "yyyy.mm"))
  filename_met = searchdir(folder_met, settings["met_prefix"] * "DATA_" * Dates.format(t, "yyyymmddHHMM"))
  met_single = matread(joinpath(folder_met, filename_met[1]))

  Sdir = met_single["sdrd"]["data"]     # "direct shortwave radiation, per inclined surface area, within topography, above canopy"
  Sdif = met_single["sdfd"]["data"]     # "diffuse shortwave radiation, per horizontal surface area, within topography, above canopy"
  LW = met_single["lwrc"]["data"]       # "longwave radiation, above topography"
  Ptot = met_single["prcs"]["data"]     # "precipitation"
  Ta = met_single["tais"]["data"]       # "air temperature"
  RH = met_single["rhus"]["data"]       # "relative humidity"
  Ua = met_single["wnss"]["data"]       # "wind speed"
  Ps = met_single["pail"]["data"]       # "local air pressure"

  thres_prec = 1.04 + 273.15
  m_prec = 0.15
  Sf = 1.0 ./ (1 .+ exp.((Ta .- thres_prec) ./ m_prec)) .* Ptot
  Rf = Ptot .- Sf

  met.Sdir[:, :] .= Sdir
  met.Sdif[:, :] .= Sdif
  met.LW[:, :] .= LW
  met.Sf[:, :] .= Sf
  met.Rf[:, :] .= Rf
  met.Ta[:, :] .= Ta
  met.RH[:, :] .= RH
  met.Ua[:, :] .= Ua
  met.Ps[:, :] .= Ps

  if fsm.SNTRAN == 1
    if haskey(met_single, "wnsd")
      Ua = met_single["wnsd"]["data"]       # "wind speed"
    else
      Ua = met_single["wnss"]["data"]       # "wind speed"
    end
    if haskey(met_single, "wndd")
      Udir = met_single["wndd"]["data"]     # "wind direction"
    else
      Udir = met_single["wnds"]["data"]     # "wind direction"
    end
    met.Ua[:, :] .= Ua
    met.Udir[:, :] .= Udir
    Udir[isnan.(Udir)] .= 0
  end

  # Update 24-hour snowfall tracking
  curr_hour = Dates.value(Hour(t)) + 1
  met.Sf24h_f64 .+= Sf
  met.Sf24h_f64 .-= met.Sf_history_f64[:, :, curr_hour]
  met.Sf_history_f64[:, :, curr_hour] = Sf
  met.Sf24h[:, :] .= met.Sf24h_f64

  # Forest-specific time-varying transmissivity data loading
  if fsm.TILE == "forest"
    folder_tvt = joinpath(settings["tvt_folder"], "YYYY." * Dates.format(t, "mm"))
    filename_tvt = searchdir(folder_tvt, "CANRAD_" * Dates.format(t, "mmddHHMM"))
    tvt_single = matread(joinpath(folder_tvt, filename_tvt[1]))
    met.Tv[:, :] .= tvt_single["stdx"]["data"]
  end

end
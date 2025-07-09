function drive(fsm, data)

  @unpack dt = fsm

  tmp = parse.(Float64, split(data))

  year = tmp[1]
  month = tmp[2]
  day = tmp[3]
  hour = tmp[4]
  Sdir[:, :] .= tmp[5]
  Sdif[:, :] .= tmp[6]
  LW[:, :] .= tmp[7]
  Sf[:, :] .= tmp[8]
  Rf[:, :] .= tmp[9]
  Ta[:, :] .= tmp[10]
  RH[:, :] .= tmp[11]
  Ua[:, :] .= tmp[12]
  Ps[:, :] .= tmp[13]
  Sf24h[:, :] .= tmp[14]

  Tc = similar(Ta)   #### hack
  es = similar(Ta)   #### hack
  # Qa = similar(Ta)   #### hack

  Ua .= max.(Ua, 0.1)

  Sf .= Sf ./ dt
  Rf .= Rf ./ dt
  Tc .= Ta .- Tm
  es .= e0 * exp.(17.5043 * Tc ./ (241.3 .+ Tc))
  Qa .= (RH ./ 100) .* eps_fsm .* es ./ Ps

  return year, month, day, hour

end


function drive!(fsm::FSM, meteo::MET{Tf,Ti}) where {Tf<:Real,Ti<:Integer}

  @unpack dt = fsm

  @unpack Tc, es, Qa, Ua, Sf, Rf, Ta, RH, Ps = meteo

  Ua .= Tf(0.7) .* Ua
  Ua .= max.(Ua, Tf(0.1))

  Sf .= Sf ./ dt
  Rf .= Rf ./ dt
  Tc .= Ta .- Tm
  es .= e0 .* exp.(Tf(17.5043) .* Tc ./ (Tf(241.3) .+ Tc))
  Qa .= (RH ./ 100) .* eps_fsm .* es ./ Ps

end

searchdir(path, key) = filter(x -> occursin(key, x), readdir(path))

function compute_psolid!(Sf, ptot, ta, thres_prec=274.19, m_prec=0.1500)

  p_corr = 1.0
  Sf .= ptot .* p_corr ./ (1 .+ exp.((ta .- thres_prec) ./ m_prec))

end

function compute_pliquid!(Rf, ptot, ta, thres_prec=274.19, m_prec=0.1500)

  p_corr = 1.0 
  Rf .= ptot .* p_corr .* exp.((ta .- thres_prec) ./ m_prec) ./ (1 .+ exp.((ta .- thres_prec) ./ m_prec))

end

function drive_grid!(meteo::MET{Tf, Ti}, fsm::FSM{Tf, Ti}, t::DateTime) where {Tf<:Real, Ti<:Integer}

  @unpack dt = fsm

  @unpack Sdir, Sdif, Sdird, LW, Sf, Rf, Ta, RH, Ua, Ps, Tc, es, Qa, Sf24h, Sf_history = meteo

  folder = joinpath("K:/DATA_ICON/OUTPUT_OSHD_0250/PROCESSED_ANALYSIS/ICON_1EFA", Dates.format(t, "yyyy.mm"))
  filename = searchdir(folder, "ICONDATA_" * Dates.format(t, "yyyymmddHHMM") * "_C1EFA_")

  met_single = matread(joinpath(folder, filename[1]))

  Sdir .= met_single["sdri"]["data"]     # "direct shortwave radiation, per inclined surface area, within topography, above canopy"
  Sdif .= met_single["sdfd"]["data"]     # "diffuse shortwave radiation, per horizontal surface area, within topography, above canopy"
  Sdird .= met_single["sdrd"]["data"]    # "direct shortwave radiation, per horizontal surface area, within topography, above canopy"
  LW .= met_single["lwrs"]["data"]       # "longwave radiation, above topography"
  Rf .= met_single["prfc"]["data"]       # "rainfall"
  Sf .= met_single["psfc"]["data"]       # "snowfall (snow + graupel)"
  Ta .= met_single["tais"]["data"]       # "air temperature"
  RH .= met_single["rhus"]["data"]       # "relative humidity"
  Ua .= met_single["wnss"]["data"]       # "wind speed"
  Ps .= met_single["pail"]["data"]       # "local air pressure"

  Ua .= Tf(0.7) .* Ua
  Ua .= max.(Ua, Tf(0.1))
  Sf .= Sf ./ dt
  Rf .= Rf ./ dt
  Tc .= Ta .- Tm
  es .= e0 .* exp.(Tf(17.5043) .* Tc ./ (Tf(241.3) .+ Tc))
  Qa .= (RH ./ 100) .* eps_fsm .* es ./ Ps

  curr_hour = Dates.value(Hour(t)) + 1 #0h -> 1; 23h -> 24
  Sf24h .+= Sf .* dt
  Sf24h .-= Sf_history[:,:,curr_hour] .* dt
  Sf_history[:,:,curr_hour] = Sf

  #TODO: missing: Tv

end

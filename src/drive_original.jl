function open_files(base_folder, fsm::FSM)

  bin_folder = joinpath(base_folder, "FSM_HS_single/bin_files")

  io = Dict()
  io["year"] = open(joinpath(bin_folder, "drive_year.bin"), "r")
  io["month"] = open(joinpath(bin_folder, "drive_month.bin"), "r")
  io["day"] = open(joinpath(bin_folder, "drive_day.bin"), "r")
  io["hour"] = open(joinpath(bin_folder, "drive_hour.bin"), "r")
  io["sdrx"] = open(joinpath(bin_folder, "drive_sdrx.bin"), "r")
  io["sdfx"] = open(joinpath(bin_folder, "drive_sdfx.bin"), "r")
  io["lwrx"] = open(joinpath(bin_folder, "drive_lwrx.bin"), "r")
  io["snfx"] = open(joinpath(bin_folder, "drive_snfx.bin"), "r")
  io["rnfx"] = open(joinpath(bin_folder, "drive_rnfx.bin"), "r")
  io["taix"] = open(joinpath(bin_folder, "drive_taix.bin"), "r")
  io["rhux"] = open(joinpath(bin_folder, "drive_rhux.bin"), "r")
  io["wnsx"] = open(joinpath(bin_folder, "drive_wnsx.bin"), "r")
  io["paix"] = open(joinpath(bin_folder, "drive_paix.bin"), "r")
  io["snfh"] = open(joinpath(bin_folder, "drive_snfh.bin"), "r")
  if (fsm.CANMOD == 1)
    io["stdx"] = open(joinpath(bin_folder, "drive_stdx.bin"), "r")
  end
  if ((fsm.ALRADT == 1) || (fsm.OSHDTN == 1))
    io["sddx"] = open(joinpath(bin_folder, "drive_sddx.bin"), "r")
  end

  return io

end


function close_files(io)
  for key in keys(io)
    close(io[key])
  end
end


function drive_binfiles!(io, meteo::MET{Tf,Ti}, fsm::FSM) where {Tf<:Real,Ti<:Integer}

  read!(io["year"], meteo.year)
  read!(io["month"], meteo.month)
  read!(io["day"], meteo.day)
  read!(io["hour"], meteo.hour)
  read!(io["sdrx"], meteo.Sdir)
  read!(io["sdfx"], meteo.Sdif)
  read!(io["lwrx"], meteo.LW)
  read!(io["snfx"], meteo.Sf)
  read!(io["rnfx"], meteo.Rf)
  read!(io["taix"], meteo.Ta)
  read!(io["rhux"], meteo.RH)
  read!(io["wnsx"], meteo.Ua)
  read!(io["paix"], meteo.Ps)
  read!(io["snfh"], meteo.Sf24h)
  if (fsm.CANMOD == 1)
    read!(io["stdx"], meteo.Tv)
  end
  if ((fsm.ALRADT == 1) || (fsm.OSHDTN == 1))
    read!(io["sddx"], meteo.Sdird)
  end

  meteo.Ua .= max.(meteo.Ua, Tf(0.1))

  for j = 1:fsm.Ny
    for i = 1:fsm.Nx

      # Convert hourly precipitation to kg/m2/s
      meteo.Sf[i, j] = meteo.Sf[i, j] / fsm.dt
      meteo.Rf[i, j] = meteo.Rf[i, j] / fsm.dt
      meteo.Tc[i, j] = meteo.Ta[i, j] - Tf(Tm)
      meteo.es[i, j] = Tf(e0) * exp(Tf(17.5043) * meteo.Tc[i, j] / (Tf(241.3) + meteo.Tc[i, j]))
      meteo.Qa[i, j] = (meteo.RH[i, j] / 100) * Tf(eps_fsm) * meteo.es[i, j] / meteo.Ps[i, j]

      # ! Safety check for Udir: in WindNinja outputs, if Ua=0, then Udir=NaN
      # ! which creates an error in SNOWTRAN3D
      # if (SNTRAN == 1) then
      #   if (isnan(Udir(i,j))) then
      #     Udir(i,j) = 0.0
      #   end if
      # end if

    end
  end

end


# searchdir(path, key) = filter(x -> occursin(key, x), readdir(path))

function drive_matfiles!(time, meteo::MET{Tf,Ti}, fsm::FSM, is_domain::BitMatrix) where {Tf<:Real,Ti<:Integer}

  folder = joinpath("K:/DATA_ICON/OUTPUT_OSHD_STAT/PROCESSED_ANALYSIS/ICON_1EFA", Dates.format(time, "yyyy.mm"))
  filename = searchdir(folder, "ICONDATA_" * Dates.format(time, "yyyymmddHHMM") * "_C1EFA_")

  met_single = matread(joinpath(folder, filename[1]))

  meteo.year .= year(time)
  meteo.month .= month(time)
  meteo.day .= day(time)
  meteo.hour .= hour(time)
  meteo.Sdir .= met_single["sdri"]["data"][:][is_domain]     # "direct shortwave radiation, per inclined surface area, within topography, above canopy"
  meteo.Sdif .= met_single["sdfd"]["data"][:][is_domain]     # "diffuse shortwave radiation, per horizontal surface area, within topography, above canopy"
  meteo.LW .= met_single["lwrs"]["data"][:][is_domain]       # "longwave radiation, above topography"
  meteo.Sf .= met_single["psfc"]["data"][:][is_domain]       # "snowfall (snow + graupel)"
  meteo.Rf .= met_single["prfc"]["data"][:][is_domain]       # "rainfall"
  meteo.Ta .= met_single["tais"]["data"][:][is_domain]       # "air temperature"
  meteo.RH .= met_single["rhus"]["data"][:][is_domain]       # "relative humidity"
  meteo.Ua .= met_single["wnss"]["data"][:][is_domain]       # "wind speed"
  meteo.Ps .= met_single["pail"]["data"][:][is_domain]       # "local air pressure"
  # if (fsm.CANMOD == 1)
  #   read!(io["stdx"], meteo.Tv)
  # end
  if ((fsm.ALRADT == 1) || (fsm.OSHDTN == 1))
    meteo.Sdird .= met_single["sdrd"]["data"][:][is_domain]    # "direct shortwave radiation, per horizontal surface area, within topography, above canopy"
  end

  meteo.Ua .= max.(meteo.Ua, Tf(0.1))
  meteo.Ua .= Tf(0.7) .* meteo.Ua

  for j = 1:fsm.Ny
    for i = 1:fsm.Nx

      # Convert hourly precipitation to kg/m2/s
      meteo.Sf[i, j] = meteo.Sf[i, j] / fsm.dt
      meteo.Rf[i, j] = meteo.Rf[i, j] / fsm.dt
      meteo.Tc[i, j] = meteo.Ta[i, j] - Tf(Tm)
      meteo.es[i, j] = Tf(e0) * exp(Tf(17.5043) * meteo.Tc[i, j] / (Tf(241.3) + meteo.Tc[i, j]))
      meteo.Qa[i, j] = (meteo.RH[i, j] / 100) * Tf(eps_fsm) * meteo.es[i, j] / meteo.Ps[i, j]

      # ! Safety check for Udir: in WindNinja outputs, if Ua=0, then Udir=NaN
      # ! which creates an error in SNOWTRAN3D
      # if (SNTRAN == 1) then
      #   if (isnan(Udir(i,j))) then
      #     Udir(i,j) = 0.0
      #   end if
      # end if

    end
  end
  
  #computation of Sf24h assuming hourly input
  curr_hour = Dates.value(Hour(time)) + 1
  meteo.Sf24h .+= meteo.Sf .* fsm.dt
  meteo.Sf24h .-= meteo.Sf_history[:,:,curr_hour] .* fsm.dt
  meteo.Sf_history[:,:,curr_hour] = meteo.Sf
  
end

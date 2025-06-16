function open_files(base_folder, fsm::FSM)

  bin_folder = joinpath(base_folder, "FSM_HS/bin_files")

  io = Dict()
  io["year"]  = open(joinpath(bin_folder, "drive_year.bin"), "r")
  io["month"] = open(joinpath(bin_folder, "drive_month.bin"), "r")
  io["day"]   = open(joinpath(bin_folder, "drive_day.bin"), "r")
  io["hour"]  = open(joinpath(bin_folder, "drive_hour.bin"), "r")
  io["sdrx"]  = open(joinpath(bin_folder, "drive_sdrx.bin"), "r")
  io["sdfx"]  = open(joinpath(bin_folder, "drive_sdfx.bin"), "r")
  io["lwrx"]  = open(joinpath(bin_folder, "drive_lwrx.bin"), "r")
  io["snfx"]  = open(joinpath(bin_folder, "drive_snfx.bin"), "r")
  io["rnfx"]  = open(joinpath(bin_folder, "drive_rnfx.bin"), "r")
  io["taix"]  = open(joinpath(bin_folder, "drive_taix.bin"), "r")
  io["rhux"]  = open(joinpath(bin_folder, "drive_rhux.bin"), "r")
  io["wnsx"]  = open(joinpath(bin_folder, "drive_wnsx.bin"), "r")
  io["paix"]  = open(joinpath(bin_folder, "drive_paix.bin"), "r")
  io["snfh"]  = open(joinpath(bin_folder, "drive_snfh.bin"), "r")
  if (fsm.CANMOD == 1)
    io["stdx"]  = open(joinpath(bin_folder, "drive_stdx.bin"), "r")
  end
  if ((fsm.ALRADT==1) || (fsm.OSHDTN == 1))
    io["sddx"]  = open(joinpath(bin_folder, "drive_sddx.bin"), "r")
  end

  return io

end


function close_files(io)
  for key in keys(io)
    close(io[key])
  end
end


function drive_original!(io, meteo::MET{Tf, Ti}, fsm::FSM) where {Tf <: Real, Ti <: Integer}

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
      meteo.Sf[i,j] = meteo.Sf[i,j]/fsm.dt
      meteo.Rf[i,j] = meteo.Rf[i,j]/fsm.dt
      meteo.Tc[i,j] = meteo.Ta[i,j] - Tf(Tm)
      meteo.es[i,j] = Tf(e0)*exp(Tf(17.5043)*meteo.Tc[i,j]/(Tf(241.3) + meteo.Tc[i,j]))
      meteo.Qa[i,j] = (meteo.RH[i,j]/100)*Tf(eps_fsm)*meteo.es[i,j]/meteo.Ps[i,j]
    
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













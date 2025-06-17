function test(base_folder, file, variable, test_vel)

  file = joinpath(base_folder, "FSM_HS/bin_files", file)

  for line in eachline(open(file))
    line = strip(line)
    if startswith(line, variable)
      parts = split(line, ' ', keepempty=false)
      if length(parts) > 1
        ref_val = parse.(Float32, parts[2:end])
        atol = isapprox(ref_val, test_vel; atol=1e-7)
        rtol = isapprox(ref_val, test_vel; rtol=1e-5)
        println("Abs tol: ", atol, " | Rel tol: ", rtol, "   ", variable)
      end
    end
  end

end


function verify_setup(base_folder, fsm)
  println("Testing setup...")
  file = "test_setup.txt"
  test(base_folder, file, "albs", fsm.albs)
  test(base_folder, file, "Ds", fsm.Ds)
  test(base_folder, file, "Nsnow", fsm.Nsnow)
  test(base_folder, file, "Qcan", fsm.Qcan)
  test(base_folder, file, "Sice", fsm.Sice)
  test(base_folder, file, "Sliq", fsm.Sliq)
  test(base_folder, file, "Sveg", fsm.Sveg)
  test(base_folder, file, "Tcan", fsm.Tcan)
  test(base_folder, file, "theta", fsm.theta)
  test(base_folder, file, "Tsnow", fsm.Tsnow)
  test(base_folder, file, "Tsoil", fsm.Tsoil)
  test(base_folder, file, "Tsrf", fsm.Tsrf)
  test(base_folder, file, "fsnow", fsm.fsnow)
  test(base_folder, file, "Tveg", fsm.Tveg)
  test(base_folder, file, "snowdepthmin", fsm.snowdepthmin)
  test(base_folder, file, "snowdepthmax", fsm.snowdepthmax)
  test(base_folder, file, "snowdepthhist", fsm.snowdepthhist)
  test(base_folder, file, "swemin", fsm.swemin)
  test(base_folder, file, "swemax", fsm.swemax)
  test(base_folder, file, "swehist", fsm.swehist)
  test(base_folder, file, "fcly", fsm.fcly)
  test(base_folder, file, "b", fsm.b)
  test(base_folder, file, "hcap_soil", fsm.hcap_soil)
  test(base_folder, file, "sathh", fsm.sathh)
  test(base_folder, file, "Vsat", fsm.Vsat)
  test(base_folder, file, "Vcrit", fsm.Vcrit)
end


function verify_drive(base_folder, meteo)
  println("Testing meteo...")
  file = "test_meteo.txt"
  test(base_folder, file, "year", meteo.year)
  test(base_folder, file, "month", meteo.month)
  test(base_folder, file, "day", meteo.day)
  test(base_folder, file, "hour", meteo.hour)
  test(base_folder, file, "Sdir", meteo.Sdir)
  test(base_folder, file, "Sdif", meteo.Sdif)
  test(base_folder, file, "LW", meteo.LW)
  test(base_folder, file, "Sf", meteo.Sf)
  test(base_folder, file, "Rf", meteo.Rf)
  test(base_folder, file, "Ta", meteo.Ta)
  test(base_folder, file, "RH", meteo.RH)
  test(base_folder, file, "Ua", meteo.Ua)
  test(base_folder, file, "Ps", meteo.Ps)
  test(base_folder, file, "Sf24h", meteo.Sf24h)
  test(base_folder, file, "Tc", meteo.Tc)
  test(base_folder, file, "es", meteo.es)
  test(base_folder, file, "Qa", meteo.Qa)
end


function verify_radiation(base_folder, fsm)
  println("Testing radiation...")
  file = "test_radiation.txt"
  test(base_folder, file, "alb", fsm.alb)
  test(base_folder, file, "asrf_out", fsm.asrf_out)
  test(base_folder, file, "Sdirt", fsm.Sdirt)
  test(base_folder, file, "Sdift", fsm.Sdift)
  test(base_folder, file, "SWveg", fsm.SWveg)
  test(base_folder, file, "SWsrf", fsm.SWsrf)
  test(base_folder, file, "SWsci", fsm.SWsci)
  test(base_folder, file, "LWt", fsm.LWt)
end


function verify_thermal(base_folder, fsm)
  println("Testing thermal...")
  file = "test_thermal.txt"
  test(base_folder, file, "Ds1 ", fsm.Ds1)
  test(base_folder, file, "gs1 ", fsm.gs1)
  test(base_folder, file, "ks1 ", fsm.ks1)
  test(base_folder, file, "Ts1 ", fsm.Ts1)
  test(base_folder, file, "Tveg0 ", fsm.Tveg0)
  test(base_folder, file, "csoil ", fsm.csoil)
  test(base_folder, file, "ksnow ", fsm.ksnow)
  test(base_folder, file, "ksoil ", fsm.ksoil )
end


function verify_sfexch(base_folder, fsm)
  println("Testing sfexch...")
  file = "test_sfexch.txt"
  test(base_folder, file, "KH ", fsm.KH)
  test(base_folder, file, "KHa ", fsm.KHa)
  test(base_folder, file, "KHg ", fsm.KHg)
  test(base_folder, file, "KHv ", fsm.KHv)
  test(base_folder, file, "KWg ", fsm.KWg)
  test(base_folder, file, "KWv ", fsm.KWv)
  test(base_folder, file, "Usc ", fsm.Usc)
end


function verify_ebalsrf(base_folder, fsm)
  println("Testing ebalsrf...")
  file = "test_ebalsrf.txt"
  test(base_folder, file, "Esrf ", fsm.Esrf)
  test(base_folder, file, "Eveg ", fsm.Eveg)
  test(base_folder, file, "G ", fsm.G)
  test(base_folder, file, "H ", fsm.H)
  test(base_folder, file, "Hsrf ", fsm.Hsrf)
  test(base_folder, file, "LE ", fsm.LE)
  test(base_folder, file, "LEsrf ", fsm.LEsrf)
  test(base_folder, file, "LWsci ", fsm.LWsci)
  test(base_folder, file, "LWveg ", fsm.LWveg)
  test(base_folder, file, "Melt ", fsm.Melt)
  test(base_folder, file, "Rnet ", fsm.Rnet)
  test(base_folder, file, "Rsrf ", fsm.Rsrf)
end


function verify_snow(base_folder, fsm)
  println("Testing snow...")
  file = "test_snow.txt"
  test(base_folder, file, "Gsoil ", fsm.Gsoil)
  test(base_folder, file, "Roff ", fsm.Roff)
  test(base_folder, file, "meltflux_out ", fsm.meltflux_out)
  test(base_folder, file, "Sbsrf ", fsm.Sbsrf)
  test(base_folder, file, "Roff_bare ", fsm.Roff_bare)
  test(base_folder, file, "Roff_snow ", fsm.Roff_snow)
  test(base_folder, file, "fsnow_thres ", fsm.fsnow_thres)
  test(base_folder, file, "unload ", fsm.unload)
  test(base_folder, file, "Tsnow ", fsm.Tsnow)
  test(base_folder, file, "Sice ", fsm.Sice)
  test(base_folder, file, "Sliq ", fsm.Sliq)
  test(base_folder, file, "Ds ", fsm.Ds)
end


function verify_soil(base_folder, fsm)
  println("Testing soil...")
  file = "test_soil.txt"
  test(base_folder, file, "Tsoil ", fsm.Tsoil)
end

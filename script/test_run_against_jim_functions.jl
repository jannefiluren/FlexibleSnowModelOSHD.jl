function test_run_against_jim_binfiles(config, fsm, meteo, ctime, station, matfiles=false)

  # Run original fortran code

  run_fsmoshd_original(config["base_folder"])

  # Determine is_domain

  landuse = matread("K:/OSHD_AUX/DATA_LUS/OSHD_LUS_STAT.mat")
  is_domain = landuse["acro"] .== station

  # Open files

  io = open_files(config["base_folder"], fsm)

  # Test setup

  if config["test_setup"]
    fields = ["albs", "Ds", "Nsnow", "Qcan", "Sice", "Sliq", "Sveg", "Tcan", "theta", "Tsnow", "Tsoil", "Tsrf", "fsnow", "Tveg", "snowdepthmin", "snowdepthmax", "snowdepthhist", "swemin", "swemax", "swehist", "fcly", "b", "hcap_soil", "sathh", "Vsat", "Vcrit"]
    verify_results(config["base_folder"], "setup", fsm, fields) && return (true, meteo)
  end

  # Loop over time

  for i = 1:24
    
    # Read drive

    dt = Dates.Hour(i-1)
    t = DateTime(year(ctime + dt), month(ctime + dt), day(ctime + dt), hour(ctime + dt), 00, 00)
    
    if matfiles
      drive_matfiles!(t, meteo, fsm, is_domain)
    else
      drive_binfiles!(io, meteo, fsm)
    end

    if config["test_drive"]
      fields = ["year", "month", "day", "hour", "Sdir", "Sdif", "LW", "Sf", "Rf", "Ta", "RH", "Ua", "Ps", "Sf24h", "Tc", "es", "Qa"]
      verify_results(config["base_folder"], "drive", meteo, fields, t) && return (true, meteo)
    end

    # Run radiation

    radiation!(fsm, meteo, t)

    if config["test_radiation"]
      fields = ["alb", "asrf_out", "Sdirt", "Sdift", "SWveg", "SWsrf", "SWsci", "LWt"]
      verify_results(config["base_folder"], "radiation", fsm, fields, t) && return (true, meteo)
    end

    # Run thermal

    thermal!(fsm)

    if config["test_thermal"]
      fields = ["Ds1", "gs1", "ks1", "Ts1", "Tveg0", "csoil", "ksnow", "ksoil"]
      verify_results(config["base_folder"], "thermal", fsm, fields, t) && return (true, meteo)
    end

    # Run sfexch and ebalsrf

    for i in 1:fsm.Nitr
      sfexch!(fsm, meteo)
      ebalsrf!(fsm, meteo)
    end

    if config["test_sfexch"]
      fields = ["KH", "KWg", "Usc"]     # ["KH", "KHa", "KHg", "KHv", "KWg", "KWv", "Usc"]
      verify_results(config["base_folder"], "sfexch", fsm, fields, t) && return (true, meteo)
    end

    if config["test_ebalsrf"]
      fields = ["Esrf", "Eveg", "G", "H", "Hsrf", "LE", "LEsrf", "LWsci", "LWveg", "Melt", "Rnet", "Rsrf"]
      verify_results(config["base_folder"], "ebalsrf", fsm, fields, t) && return (true, meteo)
    end

    # Run snow

    snow!(fsm, meteo, t)

    if config["test_snow"]
      fields = ["Gsoil",  "Roff",  "meltflux_out",  "Sbsrf",  "Roff_bare",  "Roff_snow",  "fsnow",  "unload",  "Tsnow",  "Sice",  "Sliq",  "Ds"]
      verify_results(config["base_folder"], "snow", fsm, fields, t) && return (true, meteo)
    end

    # Run soil

    soil!(fsm)

    if config["test_soil"]
      fields = ["Tsoil", "Gsoil"]
      verify_results(config["base_folder"], "soil", fsm, fields, t) && return (true, meteo)
    end

  end

  close_files(io)

  return false, meteo

end


function compile_fortran_code(base_folder)

  curr_dir = pwd()
  
  cd(joinpath(base_folder, "jim_operational/FSM_SOURCE_CODE/code"))
  run(`compil_FSM.bat`)
  cd(curr_dir)

  if !isdir(joinpath(base_folder, "FSM_HS_single/bin_files"))
    mkdir(joinpath(base_folder, "FSM_HS_single/bin_files"))
  end
  
  cp(joinpath(base_folder, "jim_operational/FSM_SOURCE_CODE/FSM2.exe"), joinpath(base_folder, "FSM_HS_single/bin_files/FSM2.exe"), force=true)

end


function run_fsmoshd_original(base_folder)

  curr_dir = pwd()

  cd(joinpath(base_folder, "FSM_HS_single/bin_files"))
  run(`FSM2.exe OPTIONS.nam`)
  cd(curr_dir)

end


function verify_results(base_folder, routine, julia_result, fields, t = nothing)
  if isnothing(t)
    println("Testing " * routine)
    file = "test_" * routine * ".txt"
  else
    println("Testing " * routine * " for hour " * @sprintf("%02d", hour(t)))
    file = "test_" * routine * "_" * @sprintf("%02d", hour(t)) * ".txt"
  end

  failed = false
  for field in fields
    if test(base_folder, file, field * " ", getfield(julia_result, Symbol(field)))
      failed = true
    end
  end
  return failed

end


function test(base_folder, file, variable, test_val)

  file = joinpath(base_folder, "FSM_HS_single/bin_files", file)

  failed = false
  open(file) do io
    for line in eachline(io)
      line = strip(line)
      if startswith(line, variable)
        parts = split(line, ' ', keepempty=false)
        if length(parts) > 1
          ref_val = parse.(Float32, parts[2:end])
          atol = isapprox(ref_val, test_val; atol=1e-4)
          rtol = isapprox(ref_val, test_val; rtol=1e-4)
          println("Abs tol: ", atol, " | Rel tol: ", rtol, "   ", variable, " (diff=", test_val - ref_val, ")")
          if atol == false && rtol == false
            failed = true
          end
        end
      end
    end
  end
  return failed

end

using FSMOSHD

# Helper functions

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

# Settings

base_folder = "D:/julia"

# Compile fsm in jim_operational

curr_dir = pwd()

cd(joinpath(base_folder, "jim_operational/FSM_SOURCE_CODE/code"))
run(`compil_FSM.bat`)
cd(curr_dir)

# Move executable to bin folder

mv(joinpath(base_folder, "jim_operational/FSM_SOURCE_CODE/FSM2.exe"), joinpath(base_folder, "FSM_HS/bin_files/FSM2.exe"), force=true)

# Run

cd(joinpath(base_folder, "FSM_HS/bin_files"))
run(`FSM2.exe OPTIONS.nam`)
cd(curr_dir)

# Create object from bin files

fsm = setup_original(Float32, Int32, joinpath(base_folder, "FSM_HS/bin_files"));
meteo = MET{Float32, Int32}(Nx=fsm.Nx, Ny=fsm.Ny)

# Test setup

if false
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

# Read drive

io = open_files(base_folder, fsm)

drive_original!(io, meteo, fsm)

if true
  println("Testing meteo...")
  file = "test_meteo.txt"
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




using ProgressMeter
using MAT
using Dates
using FSMOSHD


# tstart = DateTime(2024,9,1,6)
# tend = DateTime(2025,6,1,6)
# run_operational(tstart, tend)


function run_operational(tstart::DateTime, tend::DateTime, restart::Bool=false)

  # Settings

  base_folder = "D:/julia/FSM_julia_operational"
  
  # Read landuse data

  landuse = prepare_landuse("K:/OSHD_AUX/DATA_LUS/OSHD_LUS_0250.mat")

  # Setup model

  Nx, Ny = size(landuse["dem"]["data"])

  fsm = setup(Float32, Int32, landuse, Nx, Ny, SNFRAC=0)
  met = MET{Float32,Int32}(Nx=Nx, Ny=Ny)

  if restart
    read_states!(fsm, base_folder, tstart)
  end

  # Run model

  mkpath(joinpath(base_folder, "results"))
  mkpath(joinpath(base_folder, "states"))

  times = tstart:Hour(1):tend

  @showprogress "Running snow model..." for t in times

    # Run model for this timestep
    step!(fsm, met, t)

    if hour(t) == 5
      write_results(fsm, base_folder, t)
      write_states(fsm, base_folder, t)
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


function write_states(fsm, folder, t)

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

  file = joinpath(folder, "states", Dates.format(t + Dates.Hour(1), "yyyymmddHHMM") * "_states.mat")

  matwrite(file, states, ; compress=true)

end


function read_states!(fsm, folder, t)

  file = joinpath(folder, "states", Dates.format(t, "yyyymmddHHMM") * "_states.mat")

  states = matread(file)

  for variable in keys(states)
    state = states[variable]["data"]
    setfield!(fsm, Symbol(variable), state)
  end

end

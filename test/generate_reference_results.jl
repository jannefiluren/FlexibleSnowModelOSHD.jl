using FSMOSHD
using MAT
using Dates
using ProgressMeter


function create_full_reference_dataset()
    
    base_path = joinpath(dirname(@__FILE__), "test_data")
    reference_path = joinpath(base_path, "reference")
    
    # Configuration matrix: (tile, snfrac)
    configurations = [
        ("open", 0),
        ("open", 3), 
        ("open", 4),
        ("forest", 4)
    ]
    
    println("Generating reference results for $(length(configurations)) configurations...")
    
    for (tile, snfrac) in configurations
        println("\n" * "="^60)
        println("Generating reference for: tile=$tile, SNFRAC=$snfrac")
        println("="^60)
        
        # Run simulation for this configuration
        reference_results = run_snow_model(base_path, Tf=Float32, tile=tile, snfrac=snfrac)
        
        # Save results with configuration-specific filename
        save_reference_results(reference_results, reference_path, tile, snfrac)
    end
    
    println("\n" * "✅ All reference datasets created successfully!")
    println("Input data: $base_path")
    println("Reference results: $reference_path")
    
    return nothing
end


function run_snow_model(test_data_path::String; 
                                     Tf::Type=Float32, tile::String="open", snfrac::Int=0)
    
    println("Run simulation...")
    println("Precision: $Tf, SNFRAC: $snfrac, TILE: $tile")
    
    # Load landuse data
    landuse_path = joinpath(test_data_path, "landuse", "landuse.mat")
    if !isfile(landuse_path)
        error("Landuse file not found: $landuse_path")
    end
    
    println("Loading landuse data...")
    landuse = prepare_landuse(landuse_path)
    
    # Get domain dimensions
    Nx = size(landuse["dem"]["data"], 1)
    Ny = size(landuse["dem"]["data"], 2)
    println("Domain size: $(Nx)x$(Ny)")
    
    # Setup model
    fsm = setup(Tf, Int32, landuse, Nx, Ny, TILE=tile, SNFRAC=snfrac)
    met_curr = MET{Tf, Int32}(Nx=Nx, Ny=Ny)
    
    # Initialize 24h snowfall tracking
    Sf24h = zeros(Tf, Nx, Ny)
    Sf_history = zeros(Tf, Nx, Ny, 24)
    
    # Get list of forcing files
    meteo_path = joinpath(test_data_path, "meteo")
    tvt_path = joinpath(test_data_path, "tvt")
    
    meteo_files = filter(x -> endswith(x, ".mat"), readdir(meteo_path))
    tvt_files = filter(x -> endswith(x, ".mat"), readdir(tvt_path))
    
    sort!(meteo_files)
    sort!(tvt_files)
    
    if length(meteo_files) != length(tvt_files)
        error("Mismatch in meteorological ($(length(meteo_files))) and TVT ($(length(tvt_files))) files")
    end
    
    n_timesteps = length(meteo_files)
    println("Found $n_timesteps timesteps")
    
    # Preallocate arrays to store simulation results
    simulation_results = Dict{String, Any}()
    simulation_results["timestamps"] = Vector{String}(undef, n_timesteps)
    simulation_results["filenames"] = Vector{String}(undef, n_timesteps)
    
    state_vars = ["Tsrf", "Tsnow", "Sice", "Sliq", "fsnow", "albs", "Sveg", "Tveg", "Tcan"]
    flux_vars = ["SWsrf", "H", "LE", "G", "Melt", "Esrf", "Eveg", "Rnet"]
    diag_vars = ["Ds"]
    
    for var in [state_vars; flux_vars; diag_vars]
        simulation_results[var] = Array{Tf, 3}(undef, Nx, Ny, n_timesteps)
    end
    
    # Progress meter
    progress = Progress(n_timesteps, desc="Running reference simulation...")
    
    # Main simulation loop
    for (timestep, (meteo_file, tvt_file)) in enumerate(zip(meteo_files, tvt_files))
        
        # Load meteorological data
        met_data = matread(joinpath(meteo_path, meteo_file))
        tvt_data = matread(joinpath(tvt_path, tvt_file))

        # Exrtact time
        t = DateTime(met_data["time"], "yyyymmddHH")
        
        # Extract meteorological fields
        Sdir = met_data["sdri"]["data"]
        Sdif = met_data["sdfd"]["data"] 
        Sdird = met_data["sdrd"]["data"]
        LW = met_data["lwrs"]["data"]
        Rf = met_data["prfc"]["data"]
        Sf = met_data["psfc"]["data"]
        Ta = met_data["tais"]["data"]
        RH = met_data["rhus"]["data"]
        Ua = met_data["wnss"]["data"]
        Ps = met_data["pail"]["data"]
        
        # Extract time-varying transmissivity
        Tv = tvt_data["stdx"]["data"]
        
        # Update meteorological fields in met_curr
        met_curr.Sdir[:, :] .= Tf.(Sdir)
        met_curr.Sdif[:, :] .= Tf.(Sdif)
        met_curr.Sdird[:, :] .= Tf.(Sdird)
        met_curr.LW[:, :] .= Tf.(LW)
        met_curr.Sf[:, :] .= Tf.(Sf)
        met_curr.Rf[:, :] .= Tf.(Rf)
        met_curr.Ta[:, :] .= Tf.(Ta)
        met_curr.RH[:, :] .= Tf.(RH)
        met_curr.Ua[:, :] .= Tf.(Ua)
        met_curr.Ps[:, :] .= Tf.(Ps)
        met_curr.Tv[:, :] .= Tf.(Tv)
        
        # Update 24h snowfall tracking
        curr_hour = Dates.value(Hour(t)) + 1
        Sf24h .+= Sf
        Sf24h .-= Sf_history[:,:,curr_hour]
        Sf_history[:,:,curr_hour] = Sf
        met_curr.Sf24h[:, :] .= Sf24h
        
        # Run model for this timestep - different workflow for forest vs open tile
        drive!(fsm, met_curr)
        radiation!(fsm, met_curr, t)
        thermal!(fsm)
        
        if tile == "forest"
            # Forest workflow: use ebalfor + canopy
            for _ in 1:fsm.Nitr
                sfexch!(fsm, met_curr)
                ebalfor!(fsm, met_curr)
            end
            canopy!(fsm, met_curr)
        else
            # Open workflow: use ebalsrf
            for _ in 1:fsm.Nitr
                sfexch!(fsm, met_curr)
                ebalsrf!(fsm, met_curr)
            end
        end
        
        snow!(fsm, met_curr, t)
        soil!(fsm)
        
        # Store results
        simulation_results["timestamps"][timestep] = string(t)
        simulation_results["filenames"][timestep] = meteo_file
        
        # Store state variables
        simulation_results["Tsrf"][:, :, timestep] .= fsm.Tsrf
        simulation_results["Tsnow"][:, :, timestep] .= fsm.Tsnow[1,:,:]  # First snow layer
        simulation_results["Sice"][:, :, timestep] .= dropdims(sum(fsm.Sice, dims=1), dims=1)  # Total snow ice
        simulation_results["Sliq"][:, :, timestep] .= dropdims(sum(fsm.Sliq, dims=1), dims=1)  # Total snow liquid  
        simulation_results["fsnow"][:, :, timestep] .= fsm.fsnow
        simulation_results["albs"][:, :, timestep] .= fsm.albs
        simulation_results["Sveg"][:, :, timestep] .= fsm.Sveg
        simulation_results["Tveg"][:, :, timestep] .= fsm.Tveg
        simulation_results["Tcan"][:, :, timestep] .= fsm.Tcan
        
        # Store flux variables
        simulation_results["SWsrf"][:, :, timestep] .= fsm.SWsrf
        simulation_results["H"][:, :, timestep] .= fsm.H
        simulation_results["LE"][:, :, timestep] .= fsm.LE
        simulation_results["G"][:, :, timestep] .= fsm.G
        simulation_results["Melt"][:, :, timestep] .= fsm.Melt
        simulation_results["Esrf"][:, :, timestep] .= fsm.Esrf
        simulation_results["Eveg"][:, :, timestep] .= fsm.Eveg
        simulation_results["Rnet"][:, :, timestep] .= fsm.Rnet

        # Store diagnostic variables
        simulation_results["Ds"][:, :, timestep] .= dropdims(sum(fsm.Ds, dims=1), dims=1)  # Total snow depth  
        
        next!(progress)
    end
    
    finish!(progress)
    
    return simulation_results
end




function save_reference_results(reference_results::Dict, output_path::String, tile::String, snfrac::Int)
    
    mkpath(output_path)
    
    # Save reference results with tile and SNFRAC in filename
    reference_file = joinpath(output_path, "reference_results_$(tile)_SNFRAC_$(snfrac).mat")
    matwrite(reference_file, reference_results)
    
    # Save configuration for debugging
    config = Dict(
        "tile" => tile,
        "snfrac" => snfrac,
        "n_timesteps" => length(reference_results["timestamps"]),
        "n_gridpoints" => length(reference_results["Tsrf"][1]),
        "variables_tracked" => collect(keys(reference_results)),
        "simulation_date" => string(now())
    )
    
    config_file = joinpath(output_path, "reference_config_$(tile)_SNFRAC_$(snfrac).mat")
    matwrite(config_file, config)
    
    println("Reference results saved to: $reference_file")
    println("Configuration saved to: $config_file")
    
end

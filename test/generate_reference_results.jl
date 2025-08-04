using FSMOSHD
using MAT
using Dates
using ProgressMeter


function create_full_reference_dataset()
    
    base_path = joinpath(dirname(@__FILE__), "test_data")
    
    # Generate reference results
    println("Running reference simulation...")
    reference_results = run_snow_model(base_path)
    
    # Save reference results
    reference_path = joinpath(base_path, "reference")
    save_reference_results(reference_results, reference_path)
    
    println("Test dataset creation complete!")
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
    landuse = matread(landuse_path)
    
    # Process landuse data
    landuse = process_landuse_for_fsm(landuse)
    
    # Get domain dimensions
    Nx = size(landuse["dem"]["data"], 1)
    Ny = size(landuse["dem"]["data"], 2)
    println("Domain size: $(Nx)x$(Ny)")
    
    # Setup model
    fsm = setup_matfiles_grid(Tf, Int32, landuse, Nx, Ny, TILE=tile, SNFRAC=snfrac)
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
    
    # Initialize arrays to store simulation results
    simulation_results = Dict{String, Any}()
    simulation_results["timestamps"] = String[]
    simulation_results["filenames"] = String[]
    
    # Variables to track
    state_vars = ["Tsrf", "Tsnow", "Sice", "Sliq", "fsnow", "albs", "Sveg", "Tveg", "Tcan"]
    flux_vars = ["SWsrf", "H", "LE", "G", "Melt", "Esrf", "Eveg", "Rnet"]
    diag_vars = ["Ds"]
    
    for var in [state_vars; flux_vars; diag_vars]
        simulation_results[var] = []
    end
    
    # Progress meter
    progress = Progress(n_timesteps, desc="Running reference simulation...")
    
    # Main simulation loop
    for (meteo_file, tvt_file) in zip(meteo_files, tvt_files)
        
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
        
        # Run model for this timestep
        drive!(fsm, met_curr)
        radiation(fsm, met_curr, t)
        thermal(fsm)
        for _ in 1:fsm.Nitr
            sfexch(fsm, met_curr)
            ebalsrf(fsm, met_curr)
        end
        snow(fsm, met_curr, t)
        soil(fsm)
        
        # Store results
        push!(simulation_results["timestamps"], string(t))
        push!(simulation_results["filenames"], meteo_file)
        
        # Store state variables (flattened for easy comparison)
        push!(simulation_results["Tsrf"], vec(fsm.Tsrf))
        push!(simulation_results["Tsnow"], vec(fsm.Tsnow[1,:,:]))  # First snow layer
        push!(simulation_results["Sice"], vec(sum(fsm.Sice, dims=1)))  # Total snow ice
        push!(simulation_results["Sliq"], vec(sum(fsm.Sliq, dims=1)))  # Total snow liquid  
        push!(simulation_results["fsnow"], vec(fsm.fsnow))
        push!(simulation_results["albs"], vec(fsm.albs))
        push!(simulation_results["Sveg"], vec(fsm.Sveg))
        push!(simulation_results["Tveg"], vec(fsm.Tveg))
        push!(simulation_results["Tcan"], vec(fsm.Tcan))
        
        # Store flux variables
        push!(simulation_results["SWsrf"], vec(fsm.SWsrf))
        push!(simulation_results["H"], vec(fsm.H))
        push!(simulation_results["LE"], vec(fsm.LE))
        push!(simulation_results["G"], vec(fsm.G))
        push!(simulation_results["Melt"], vec(fsm.Melt))
        push!(simulation_results["Esrf"], vec(fsm.Esrf))
        push!(simulation_results["Eveg"], vec(fsm.Eveg))
        push!(simulation_results["Rnet"], vec(fsm.Rnet))

        # Store diagnostic variables
        push!(simulation_results["Ds"], vec(sum(fsm.Ds, dims=1)))  # Total snow depth  
        
        next!(progress)
    end
    
    finish!(progress)
    
    return simulation_results
end


function process_landuse_for_fsm(landuse::Dict)
    
    processed_landuse = copy(landuse)
    
    # Add required fields that are computed in prepare_landuse_grid
    dem_size = size(processed_landuse["dem"]["data"])
    
    # processed_landuse["is_domain"] = ones(Bool, dem_size)
    processed_landuse["dhdxdy"] = processed_landuse["dhdxdy"]["data"]   # TODO change in subsequent code
    processed_landuse["sddem"] = processed_landuse["sd"]["data"]   # TODO change in subsequent code
    processed_landuse["Ld"] = 250 * ones(dem_size)
    processed_landuse["slopemu"] = sqrt.((processed_landuse["dhdxdy"] ./ 2))
    processed_landuse["xi"] = (sqrt(2) * processed_landuse["sddem"]) ./ processed_landuse["slopemu"]
    processed_landuse["x"] = ones(dem_size)
    processed_landuse["y"] = ones(dem_size)
    processed_landuse["prec_multi"]["data"] .= 1
    
    return processed_landuse
end


function save_reference_results(reference_results::Dict, output_path::String)
    
    mkpath(output_path)
    
    # Save reference results
    reference_file = joinpath(output_path, "reference_results.mat")
    matwrite(reference_file, reference_results)
    
    println("Reference results saved to: $reference_file")
    
end

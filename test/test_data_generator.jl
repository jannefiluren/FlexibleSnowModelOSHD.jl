using Dates
using MAT
using FSMOSHD


function create_test_dataset(;
    start_date::Date=Date(2024, 10, 1),
    end_date::Date=Date(2025, 6, 1),
    domain_size::Int=20)

    base_path = joinpath(dirname(@__FILE__), "test_data")

    println("Creating test dataset from real data...")
    println("Period: $start_date to $end_date")
    println("Domain size: $(domain_size)x$(domain_size)")

    landuse = load_landuse_data(domain_size=domain_size)

    meteo_data, tvt_data = load_forcing_data(start_date, end_date, landuse)

    save_test_data(landuse, meteo_data, tvt_data, base_path)

    println("Test dataset creation complete!")
end


function extract_test_domain_from_landuse(landuse_data::Dict, center_x::Int, center_y::Int, domain_size::Int)

    # Calculate domain bounds
    half_size = domain_size ÷ 2
    x_start = center_x - half_size
    x_end = center_x + half_size - 1
    y_start = center_y - half_size
    y_end = center_y + half_size - 1

    # Ensure bounds are within the full dataset
    full_nx, full_ny = size(landuse_data["dem"]["data"])
    x_start = max(1, x_start)
    x_end = min(full_nx, x_end)
    y_start = max(1, y_start)
    y_end = min(full_ny, y_end)

    println("Extracting domain from ($x_start:$x_end, $y_start:$y_end)")

    # Extract all relevant landuse fields
    data_fields = ["dem", "skyvf", "fveg", "hcan", "lai", "vfhp", "fves", "prec_multi",
        "forest", "slope", "dhdxdy", "sd"]

    test_landuse = Dict()
    for field in data_fields
        if haskey(landuse_data, field) && haskey(landuse_data[field], "data")
            test_landuse[field] = Dict("data" => landuse_data[field]["data"][x_start:x_end, y_start:y_end])
        end
    end

    # Add field for glacier tile
    test_landuse["glacier"] = Dict("data" => ones(domain_size, domain_size))

    # Add domain definition
    test_landuse["is_domain"] = Dict("data" => ones(Bool, size(test_landuse["dem"]["data"])))

    # Add indicies of domain bounds
    test_landuse["x_start"] = x_start
    test_landuse["x_end"] = x_end
    test_landuse["y_start"] = y_start
    test_landuse["y_end"] = y_end

    return test_landuse
end


function load_landuse_data(; domain_size::Int=20, center_x::Union{Int,Nothing}=nothing, center_y::Union{Int,Nothing}=nothing)

    println("Loading landuse data...")

    landuse_file = "K:/OSHD_AUX/DATA_LUS/OSHD_LUS_0250.mat"

    if !isfile(landuse_file)
        error("Landuse file not found: $landuse_file. Please ensure the data path is correct.")
    end

    full_landuse = matread(landuse_file)

    # Auto-select center if not provided (choose area with interesting topography)
    full_nx, full_ny = size(full_landuse["dem"]["data"])
    if center_x === nothing || center_y === nothing
        center_x = full_nx ÷ 2
        center_y = full_ny ÷ 2
        println("Auto-selected domain center: ($center_x, $center_y)")
    end

    # Extract test domain
    test_landuse = extract_test_domain_from_landuse(full_landuse, center_x, center_y, domain_size)

    println("Extracted $(domain_size)x$(domain_size) domain from full landuse dataset")
    println("Elevation range: $(round(minimum(test_landuse["dem"]["data"]), digits=1)) - $(round(maximum(test_landuse["dem"]["data"]), digits=1)) m")

    return test_landuse
end


function load_forcing_data(start_date::Date, end_date::Date, landuse::Dict;
    base_path_icon::String="K:/DATA_ICON/OUTPUT_OSHD_0250/PROCESSED_ANALYSIS/ICON_1EFA",
    base_path_canrad::String="K:/OSHD_AUX/DATA_CANRAD/OUTPUT_OSHD_0250/CR_2410_static")

    # Extract domain
    domain_x_range = landuse["x_start"]:landuse["x_end"]
    domain_y_range = landuse["y_start"]:landuse["y_end"]

    # Create hourly datetime range
    datetimes = DateTime(start_date):Hour(1):DateTime(end_date)
    
    
    println("Loading meteorological and time-varying transmissivity data")
    println("Period: $start_date to $end_date")
    println("Total timesteps: $(length(datetimes))")
    
    meteo_data, tvt_data = [], []
    for t in datetimes

        # Load meteorological data
        folder_icon = joinpath(base_path_icon, Dates.format(t, "yyyy.mm"))
        filename_icon = searchdir(folder_icon, "ICONDATA_" * Dates.format(t, "yyyymmddHHMM") * "_C1EFA_")
        met_raw = matread(joinpath(folder_icon, filename_icon[1]))

        met_single = Dict{String, Any}()
        met_single["sdri"] = Dict("data" => met_raw["sdri"]["data"][domain_x_range, domain_y_range])
        met_single["sdfd"] = Dict("data" => met_raw["sdfd"]["data"][domain_x_range, domain_y_range])
        met_single["sdrd"] = Dict("data" => met_raw["sdrd"]["data"][domain_x_range, domain_y_range])
        met_single["lwrs"] = Dict("data" => met_raw["lwrs"]["data"][domain_x_range, domain_y_range])
        met_single["prfc"] = Dict("data" => met_raw["prfc"]["data"][domain_x_range, domain_y_range])
        met_single["psfc"] = Dict("data" => met_raw["psfc"]["data"][domain_x_range, domain_y_range])
        met_single["tais"] = Dict("data" => met_raw["tais"]["data"][domain_x_range, domain_y_range])
        met_single["rhus"] = Dict("data" => met_raw["rhus"]["data"][domain_x_range, domain_y_range])
        met_single["wnss"] = Dict("data" => met_raw["wnss"]["data"][domain_x_range, domain_y_range])
        met_single["pail"] = Dict("data" => met_raw["pail"]["data"][domain_x_range, domain_y_range])
        met_single["filename"] = filename_icon[1]
        met_single["time"] = filename_icon[1][10:19]

        push!(meteo_data, met_single)

        # Load time-varying transmissivity data
        folder_tvt = joinpath(base_path_canrad, "YYYY." * Dates.format(t, "mm"))
        filename_tvt = searchdir(folder_tvt, "CANRAD_" * Dates.format(t, "mmddHHMM"))

        tvt_raw = matread(joinpath(folder_tvt, filename_tvt[1]))
        
        tvt_single = Dict{String, Any}()
        tvt_single["stdx"] = Dict{String, Any}()
        tvt_single["stdx"]["data"] = tvt_raw["stdx"]["data"][domain_x_range, domain_y_range]
        tvt_single["stdx"]["unit"] = "-"
        tvt_single["filename"] = filename_tvt[1]

        push!(tvt_data, tvt_single)

    end

    println("Successfully loaded $(length(meteo_data)) files")
    return meteo_data, tvt_data
end


function save_test_data(landuse, meteo_data, tvt_data, base_path::String)
    
    println("Saving datasets...")

    # Save landuse data
    mkpath(joinpath(base_path, "landuse"))
    matwrite(joinpath(base_path, "landuse", "landuse.mat"), landuse)

    # Save meteorology data
    mkpath(joinpath(base_path, "meteo"))
    for met in meteo_data
        filename = met["filename"]
        matwrite(joinpath(base_path, "meteo", filename), met)
    end

    # Save transmissivity data
    mkpath(joinpath(base_path, "tvt"))
    for tvt in tvt_data
        filename = tvt["filename"]
        matwrite(joinpath(base_path, "tvt", filename), tvt)
    end

    println("Test data saved to: $base_path")
end
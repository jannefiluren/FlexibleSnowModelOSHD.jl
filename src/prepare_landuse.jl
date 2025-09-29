"""
    prepare_landuse(filename::String)

Load landuse data from file and fill missing fields with default values.

# Arguments
- `filename::String`: Path to landuse .mat file

# Returns
- `Dict`: Processed landuse dictionary with all required fields and computed derivatives
"""
function prepare_landuse(filename::String)
    
    # Load landuse data
    landuse = matread(filename)

    # Handle old landuse format
    if haskey(landuse, "landuse")
        landuse = landuse["landuse"]
    end
    
    # Get domain dimensions
    if !isa(landuse["dem"], Dict)
        dem_size = size(landuse["dem"])
    else
        dem_size = size(landuse["dem"]["data"])
    end
    
    # Ensure all fields have the ["data"] format or fill missing fields with defaults
    _ensure_data_field!(landuse, "dem", ones(Bool, dem_size))
    _ensure_data_field!(landuse, "is_domain", ones(Bool, dem_size))
    _ensure_data_field!(landuse, "dhdxdy", ones(dem_size))
    _ensure_data_field!(landuse, "sd", ones(dem_size))
    _ensure_data_field!(landuse, "x", ones(dem_size))
    _ensure_data_field!(landuse, "y", ones(dem_size))
    _ensure_data_field!(landuse, "prec_multi", ones(dem_size))
    _ensure_data_field!(landuse, "skyvf", ones(dem_size))

    # Handle the special case of grid cell size (not available for station simulations)
    if haskey(landuse, "cellsize")
        _ensure_data_field!(landuse, "Ld", landuse["cellsize"] * ones(dem_size))
    else
        _ensure_data_field!(landuse, "Ld", zeros(dem_size))
    end
    
    # Handle optional forest fields
    _ensure_data_field!(landuse, "fveg", nothing)
    _ensure_data_field!(landuse, "hcan", nothing)
    _ensure_data_field!(landuse, "lai", nothing)
    _ensure_data_field!(landuse, "vfhp", nothing)
    _ensure_data_field!(landuse, "fves", nothing)
    _ensure_data_field!(landuse, "forest", nothing)
    
    # Compute derived fields using the data arrays
    landuse["slopemu"] = Dict("data" => sqrt.((landuse["dhdxdy"]["data"] ./ 2)))
    landuse["xi"] = Dict("data" => (sqrt(2) * landuse["sd"]["data"]) ./ landuse["slopemu"]["data"])
    
    return landuse
end

"""
    crop_landuse_to_domain(landuse::Dict)

Crop landuse data to only include points where is_domain is true.
This converts 2D arrays to 1D vectors containing only the selected points.

# Arguments
- `landuse::Dict`: Landuse dictionary with ["variable"]["data"] format

# Returns
- `Dict`: Cropped landuse dictionary where 2D arrays become 1D vectors
"""
function crop_landuse_to_domain(landuse::Dict)
    
    # Get the domain mask
    is_domain = landuse["is_domain"]["data"]
    
    # Create cropped landuse dictionary
    cropped_landuse = deepcopy(landuse)
    
    # List of variables that should be cropped (2D spatial arrays)
    spatial_vars = [
        "dem", "skyvf", "x", "y", "dhdxdy", "sd", "Ld", "slopemu", "xi",
        "fveg", "hcan", "lai", "vfhp", "fves", "prec_multi", "forest"
    ]
    
    # Crop each spatial variable
    for var in spatial_vars
        if haskey(cropped_landuse, var) && haskey(cropped_landuse[var], "data")
            # Extract only the domain points and reshape to column vector
            cropped_data = cropped_landuse[var]["data"][is_domain]
            cropped_landuse[var]["data"] = reshape(cropped_data, :, 1)
        end
    end
    
    # Update is_domain to be a 2D array of all trues
    cropped_landuse["is_domain"]["data"] = ones(Bool, sum(is_domain), 1)
    
    return cropped_landuse
end

"""
Helper function to ensure a landuse field has the ["data"] format
"""
function _ensure_data_field!(landuse::Dict, field::String, default_value)
    if !haskey(landuse, field)
        # Field doesn't exist - create it with default if provided
        if default_value !== nothing
            landuse[field] = Dict("data" => default_value)
        end
    elseif !isa(landuse[field], Dict)
        # Field exists but is not a Dict - wrap it in ["data"] format
        landuse[field] = Dict("data" => landuse[field])
    elseif !haskey(landuse[field], "data")
        # Field is a Dict but missing "data" key - add default if provided
        if default_value !== nothing
            landuse[field]["data"] = default_value
        end
    end
    # If field exists and already has proper ["data"] format, do nothing
end

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
    
    # Get domain dimensions
    dem_size = size(landuse["dem"]["data"])
    
    # Ensure all fields have the ["data"] format or fill missing fields with defaults
    _ensure_data_field!(landuse, "is_domain", ones(Bool, dem_size))
    _ensure_data_field!(landuse, "dhdxdy", ones(dem_size))
    _ensure_data_field!(landuse, "sd", ones(dem_size))
    _ensure_data_field!(landuse, "Ld", 250 * ones(dem_size))
    _ensure_data_field!(landuse, "x", ones(dem_size))
    _ensure_data_field!(landuse, "y", ones(dem_size))
    _ensure_data_field!(landuse, "prec_multi", ones(dem_size))
    _ensure_data_field!(landuse, "skyvf", ones(dem_size))
    
    # Handle optional forest fields
    _ensure_data_field!(landuse, "fveg", nothing)
    _ensure_data_field!(landuse, "hcan", nothing)
    _ensure_data_field!(landuse, "lai", nothing)
    _ensure_data_field!(landuse, "vfhp", nothing)
    _ensure_data_field!(landuse, "fves", nothing)
    _ensure_data_field!(landuse, "forest", nothing)
    
    # Set precipitation multiplier to one    TODO handle this later, in particular for the forest case...
    landuse["prec_multi"]["data"] .= 1
    
    # Compute derived fields using the data arrays
    landuse["slopemu"] = Dict("data" => sqrt.((landuse["dhdxdy"]["data"] ./ 2)))
    landuse["xi"] = Dict("data" => (sqrt(2) * landuse["sd"]["data"]) ./ landuse["slopemu"]["data"])
    
    return landuse
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

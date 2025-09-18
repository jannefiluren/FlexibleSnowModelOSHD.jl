using MAT
using Dates

# Variables available for saving to file

const AVAILABLE_OUTPUT_VARS = Dict(
    "Roff" => Dict(
        "longname" => "total runoff",
        "acronymn" => "rotc",
        "unit" => "mm/tstep",
        "type" => "sum"
        ),
    "Roff_snow" => Dict(
        "longname" => "runoff from snow",
        "acronymn" => "romc",
        "unit" => "mm/tstep",
        "type" => "sum"
        ),
    "Sbsrf" => Dict(
        "longname" => "sublimation from snow",
        "acronymn" => "sbsc",
        "unit" => "mm/tstep",
        "type" => "sum"
        ),
    "fsnow" => Dict(
        "longname" => "snow covered fraction",
        "acronymn" => "scfe",
        "unit" => "-",
        "type" => "instant"
        ),
    "snowdepth" => Dict(
        "longname" => "snow depth",
        "acronymn" => "hsnt",
        "unit" => "m",
        "type" => "instant"
        ),
    "SWE" => Dict(
        "longname" => "snow water equivalent",
        "acronymn" => "swet",
        "unit" => "mm",
        "type" => "instant"
        ),
    "Tsrf" => Dict(
        "longname" => "surface temperature",
        "acronymn" => "tsfe",
        "unit" => "K",
        "type" => "instant"
        ),
    "Sdirt" => Dict(
        "longname" => "incoming direct shortwave radiation",
        "acronymn" => "swtb",
        "unit" => "W/m^2",
        "type" => "mean"
        ),
    "Sdift" => Dict(
        "longname" => "incoming diffuse shortwave radiation",
        "acronymn" => "swtd",
        "unit" => "W/m^2",
        "type" => "mean"
        ),
    "LWt" => Dict(
        "longname" => "incoming longwave radiation",
        "acronymn" => "lwtr",
        "unit" => "W/m^2",
        "type" => "mean"
        ),
    "Sliq_out" => Dict(
        "longname" => "liquid water content",
        "acronymn" => "slqt",
        "unit" => "mm",
        "type" => "instant"
        ),
    "Nsnow" => Dict(
        "longname" => "number of snow layers",
        "acronymn" => "nsne",
        "unit" => "-",
        "type" => "instant"
        ),
    "albs" => Dict(
        "longname" => "snow albedo",
        "acronymn" => "alse",
        "unit" => "-",
        "type" => "instant"
        ),
    "H" => Dict(
        "longname" => "sensible heat flux",
        "acronymn" => "sehe",
        "unit" => "W/m^2",
        "type" => "mean"
        ),
    "LE" => Dict(
        "longname" => "latent heat flux",
        "acronymn" => "lahe",
        "unit" => "W/m^2",
        "type" => "mean"
        ),
    "Rnet" => Dict(
        "longname" => "net radiation",
        "acronymn" => "rnet",
        "unit" => "W/m^2",
        "type" => "mean"
        ),
    "asrf_out" => Dict(
        "longname" => "surface albedo",
        "acronymn" => "asrf",
        "unit" => "-",
        "type" => "instant"
        ),
    "Tsnow" => Dict(
        "longname" => "snow layer temperatures",
        "acronymn" => "tsnl",
        "unit" => "K",
        "type" => "instant"
        ),
    "Tsoil" => Dict(
        "longname" => "soil layer temperatures",
        "acronymn" => "tsll",
        "unit" => "K",
        "type" => "instant"
        ),
    "Tcan" => Dict(
        "longname" => "canopy air space temperature",
        "acronymn" => "tcan",
        "unit" => "K",
        "type" => "instant"
        ),
    "Tveg" => Dict(
        "longname" => "vegetation temperature",
        "acronymn" => "tveg",
        "unit" => "K",
        "type" => "instant"
        ),
    "Qcan" => Dict(
        "longname" => "canopy air space humidity",
        "acronymn" => "qcan",
        "unit" => "-",    # TODO what is the unit of this?
        "type" => "instant"
        ),
    "Sveg" => Dict(
        "longname" => "snow intercepted in the canopy",
        "acronymn" => "sveg",
        "unit" => "mm",
        "type" => "instant"
        ),
    "SWsrf" => Dict(
        "longname" => "net shortwave radiation absorbed by the surface",
        "acronymn" => "swsr",
        "unit" => "W/m^2",
        "type" => "mean"
        ),
    "SWveg" => Dict(
        "longname" => "net shortwave radiation absorbed by the vegetation",
        "acronymn" => "swve",
        "unit" => "W/m^2",
        "type" => "mean"
        ),
    "unload" => Dict(
        "longname" => "snow mass unloaded from canopy",
        "acronymn" => "unld",
        "unit" => "mm/tstep",
        "type" => "mean"
        ),
    "intcpt" => Dict(
        "longname" => "snow mass intercepted in the canopy",
        "acronymn" => "intc",
        "unit" => "mm/tstep",
        "type" => "mean"
        ),
    "Esrf" => Dict(
        "longname" => "moisture flux from the surface",
        "acronymn" => "esrf",
        "unit" => "kg/m^2/s",
        "type" => "mean"
        ),
    "Eveg" => Dict(
        "longname" => "moisture flux from vegetation",
        "acronymn" => "eveg",
        "unit" => "kg/m^2/s",
        "type" => "mean"
        ),
    "Gsoil" => Dict(
        "longname" => "heat flux into soil",
        "acronymn" => "ghsl",
        "unit" => "W/m^2",
        "type" => "mean"
        ),
    "Hsrf" => Dict(
        "longname" => "sensible heat flux from the surface",
        "acronymn" => "hesr",
        "unit" => "W/m^2",
        "type" => "mean"
        ),
    "LEsrf" => Dict(
        "longname" => "latent heat flux from the surface",
        "acronymn" => "lesr",
        "unit" => "W/m^2",
        "type" => "mean"
        ),
    "Rsrf" => Dict(
        "longname" => "net radiation at surface",
        "acronymn" => "rnsr",
        "unit" => "W/m^2",
        "type" => "mean"
        ),
    "Melt" => Dict(
        "longname" => "surface melt rate",
        "acronymn" => "emlt",
        "unit" => "kg/m^2/s",
        "type" => "mean"
        ),
    "SWsci" => Dict(
        "longname" => "subcanopy incoming shortwave radiation",
        "acronymn" => "swsc",
        "unit" => "W/m^2",
        "type" => "mean"
        ),
    "LWsci" => Dict(
        "longname" => "subcanopy incoming longwave radiation",
        "acronymn" => "lwsc",
        "unit" => "W/m^2",
        "type" => "mean"
        ),
    "LWveg" => Dict(
        "longname" => "net longwave radiation absorbed by vegetation",
        "acronymn" => "lwve",
        "unit" => "W/m^2",
        "type" => "mean"
        ),
    "Usc" => Dict(
        "longname" => "wind speed in canopy layer",
        "acronymn" => "uasc",
        "unit" => "m/s",
        "type" => "mean"
        ),
    "Sbveg" => Dict(
        "longname" => "sublimation from vegetation",
        "acronymn" => "sbve",
        "unit" => "mm/tstep",
        "type" => "sum"
        ),
    "KH" => Dict(
        "longname" => "Eddy diffusivity for heat to the atmosphere",
        "acronymn" => "khag",
        "unit" => "m/s",
        "type" => "mean"
        ),
    "KHa" => Dict(
        "longname" => "Eddy diffusivity for heat from the canopy air space",
        "acronymn" => "khac",
        "unit" => "m/s",
        "type" => "mean"
        ),
    "KHg" => Dict(
        "longname" => "Eddy diffusivity for heat from the ground",
        "acronymn" => "khgr",
        "unit" => "m/s",
        "type" => "mean"
        ),
    "KHv" => Dict(
        "longname" => "Eddy diffusivity for heat from vegetation",
        "acronymn" => "khve",
        "unit" => "m/s",
        "type" => "mean"
        ),
    "KWg" => Dict(
        "longname" => "Eddy diffusivity for water from the ground",
        "acronymn" => "kwgr",
        "unit" => "m/s",
        "type" => "mean"
        ),
    "KWv" => Dict(
        "longname" => "Eddy diffusivity for water from vegetation",
        "acronymn" => "kwve",
        "unit" => "m/s",
        "type" => "mean"
        ),
    "snowdepthmin" => Dict(
        "longname" => "minimum snow depth at time step of swemin",
        "acronymn" => "hsmn",
        "unit" => "m",
        "type" => "instant"
        ),
    "snowdepthmax" => Dict(
        "longname" => "maximum snow depth at time stemp of swemax",
        "acronymn" => "hsmx",
        "unit" => "m",
        "type" => "instant"
        ),
    "snowdepthhist" => Dict(
        "longname" => "history of snow depth during last 14 days (most recent entries first)",
        "acronymn" => "hshs",
        "unit" => "m",
        "type" => "instant"
        ),
    "swemin" => Dict(
        "longname" => "minimum swe during the season",
        "acronymn" => "swmn",
        "unit" => "mm",
        "type" => "instant"
        ),
    "swemax" => Dict(
        "longname" => "maximum swe during the season",
        "acronymn" => "swmx",
        "unit" => "mm",
        "type" => "instant"
        ),
    "swehist" => Dict(
        "longname" => "history of SWE during last 14 days",
        "acronymn" => "swhs",
        "unit" => "mm",
        "type" => "instant"
        )
)


function save_model_output!(output_dict::Dict, fsm::FSM, t::DateTime, filepath::String; compress::Bool=true)
    
    # Nothing to save
    isempty(output_dict) && return

    # Assign data to dict fields      # TODO add snowdepth, SWE and Sliq_out as diagnostics to fsm struct?
    for var_name in keys(output_dict)
        if isa(output_dict[var_name], Dict)
            if var_name == "snowdepth"
                output_dict[var_name]["data"] = dropdims(sum(fsm.Ds, dims=1), dims=1) .* fsm.fsnow
            elseif var_name == "SWE"
                output_dict[var_name]["data"] = dropdims(sum(fsm.Sice .+ fsm.Sliq, dims=1), dims=1)
            elseif var_name == "Sliq_out"
                output_dict[var_name]["data"] = dropdims(sum(fsm.Sliq, dims=1), dims=1)
            elseif startswith(output_dict[var_name]["type"], "instant") 
                output_dict[var_name]["data"] = getfield(fsm, Symbol(var_name))
            elseif startswith(output_dict[var_name]["type"], "mean")
                output_dict[var_name]["data"] ./= output_dict["hours_in_period"]
            end
        end
    end

    # Add metadata
    output_dict["tend"] = Dates.format(t + Second(fsm.dt), "yyyy-mm-ddTHH:MM:SS")

    # Save to file
    try
        matwrite(filepath, output_dict; compress=compress)
    catch e
        error("Failed to write output file $filepath: $e")
    end

    # Reset cumulative variables to zero
    for var_name in keys(output_dict)
        if isa(output_dict[var_name], Dict)
            if startswith(output_dict[var_name]["type"], "sum") || startswith(output_dict[var_name]["type"], "mean")
                output_dict[var_name]["data"] .= 0.0
            end
        end
    end
    output_dict["hours_in_period"] = 0

    # Set start time
    output_dict["tstart"] = Dates.format(t + Second(fsm.dt), "yyyy-mm-ddTHH:MM:SS")

end


function allocate_output_dict(fsm::FSM, t::DateTime, output_vars::Vector{String})

    invalid_vars = setdiff(output_vars, keys(AVAILABLE_OUTPUT_VARS))
    if !isempty(invalid_vars)
        error("Invalid output variables: $(join(invalid_vars, ", ")). Available variables: $(join(keys(AVAILABLE_OUTPUT_VARS), ", "))")
    end

    output_dict = Dict()

    for output_var in output_vars

        size_var = (fsm.Nx, fsm.Ny)
        if Symbol(output_var) in fieldnames(typeof(fsm))
            size_var = size(getfield(fsm, Symbol(output_var)))
        end
        
        type_str = AVAILABLE_OUTPUT_VARS[output_var]["type"] * " from tstart to tend"
        if AVAILABLE_OUTPUT_VARS[output_var]["type"] == "instant"
            type_str = "instant value at tend"
        end

        output_dict[output_var] = Dict(
            "data" => fill(0.0, size_var),
            "shortname" => output_var,
            "longname" => AVAILABLE_OUTPUT_VARS[output_var]["longname"],
            "acronymn" => AVAILABLE_OUTPUT_VARS[output_var]["acronymn"],
            "unit" => AVAILABLE_OUTPUT_VARS[output_var]["unit"],
            "type" => type_str
            )

    end

    if !isempty(output_dict)
        output_dict["tstart"] = Dates.format(t, "yyyy-mm-ddTHH:MM:SS")
        output_dict["hours_in_period"] = 0
    end

    return output_dict

end


function cumulate!(output_dict::Dict, fsm::FSM)

    isempty(output_dict) && return

    for var_name in keys(output_dict)
        if isa(output_dict[var_name], Dict)
            if startswith(output_dict[var_name]["type"], "sum") || startswith(output_dict[var_name]["type"], "mean")
                output_dict[var_name]["data"] .+= getfield(fsm, Symbol(var_name))
            end
        end
    end

    output_dict["hours_in_period"] += 1

end

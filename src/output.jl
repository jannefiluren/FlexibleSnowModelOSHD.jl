@enum AggregationType INSTANT SUM MEAN

# Variables available for saving to file

const AVAILABLE_OUTPUT_VARS = Dict(
    "Roff" => Dict(
        "longname" => "total runoff",
        "acronymn" => "rotc",
        "unit" => "mm/tstep",
        "type" => SUM
    ),
    "meltflux_out" => Dict(
        "longname" => "runoff from snow",
        "acronymn" => "romc",
        "unit" => "mm/tstep",
        "type" => SUM
    ),
    "Sbsrf" => Dict(
        "longname" => "sublimation from snow",
        "acronymn" => "sbsc",
        "unit" => "mm/tstep",
        "type" => SUM
    ),
    "fsnow" => Dict(
        "longname" => "snow covered fraction",
        "acronymn" => "scfe",
        "unit" => "-",
        "type" => INSTANT
    ),
    "snowdepth" => Dict(
        "longname" => "snow depth",
        "acronymn" => "hsnt",
        "unit" => "m",
        "type" => INSTANT
    ),
    "SWE" => Dict(
        "longname" => "snow water equivalent",
        "acronymn" => "swet",
        "unit" => "mm",
        "type" => INSTANT
    ),
    "Tsrf" => Dict(
        "longname" => "surface temperature",
        "acronymn" => "tsfe",
        "unit" => "K",
        "type" => INSTANT
    ),
    "Sdirt" => Dict(
        "longname" => "incoming direct shortwave radiation",
        "acronymn" => "swtb",
        "unit" => "W/m^2",
        "type" => MEAN
    ),
    "Sdift" => Dict(
        "longname" => "incoming diffuse shortwave radiation",
        "acronymn" => "swtd",
        "unit" => "W/m^2",
        "type" => MEAN
    ),
    "LWt" => Dict(
        "longname" => "incoming longwave radiation",
        "acronymn" => "lwtr",
        "unit" => "W/m^2",
        "type" => MEAN
    ),
    "Sliq_out" => Dict(
        "longname" => "liquid water content",
        "acronymn" => "slqt",
        "unit" => "mm",
        "type" => INSTANT
    ),
    "Nsnow" => Dict(
        "longname" => "number of snow layers",
        "acronymn" => "nsne",
        "unit" => "-",
        "type" => INSTANT
    ),
    "albs" => Dict(
        "longname" => "snow albedo",
        "acronymn" => "alse",
        "unit" => "-",
        "type" => INSTANT
    ),
    "H" => Dict(
        "longname" => "sensible heat flux",
        "acronymn" => "sehe",
        "unit" => "W/m^2",
        "type" => MEAN
    ),
    "LE" => Dict(
        "longname" => "latent heat flux",
        "acronymn" => "lahe",
        "unit" => "W/m^2",
        "type" => MEAN
    ),
    "Rnet" => Dict(
        "longname" => "net radiation",
        "acronymn" => "rnet",
        "unit" => "W/m^2",
        "type" => MEAN
    ),
    "asrf_out" => Dict(
        "longname" => "surface albedo",
        "acronymn" => "asrf",
        "unit" => "-",
        "type" => INSTANT
    ),
    "Tsnow" => Dict(
        "longname" => "snow layer temperatures",
        "acronymn" => "tsnl",
        "unit" => "K",
        "type" => INSTANT
    ),
    "Tsoil" => Dict(
        "longname" => "soil layer temperatures",
        "acronymn" => "tsll",
        "unit" => "K",
        "type" => INSTANT
    ),
    "Tcan" => Dict(
        "longname" => "canopy air space temperature",
        "acronymn" => "tcan",
        "unit" => "K",
        "type" => INSTANT
    ),
    "Tveg" => Dict(
        "longname" => "vegetation temperature",
        "acronymn" => "tveg",
        "unit" => "K",
        "type" => INSTANT
    ),
    "Qcan" => Dict(
        "longname" => "canopy air space humidity",
        "acronymn" => "qcan",
        "unit" => "-",    # TODO what is the unit of this?
        "type" => INSTANT
    ),
    "Sveg" => Dict(
        "longname" => "snow intercepted in the canopy",
        "acronymn" => "sveg",
        "unit" => "mm",
        "type" => INSTANT
    ),
    "SWsrf" => Dict(
        "longname" => "net shortwave radiation absorbed by the surface",
        "acronymn" => "swsr",
        "unit" => "W/m^2",
        "type" => MEAN
    ),
    "SWveg" => Dict(
        "longname" => "net shortwave radiation absorbed by the vegetation",
        "acronymn" => "swve",
        "unit" => "W/m^2",
        "type" => MEAN
    ),
    "unload" => Dict(
        "longname" => "snow mass unloaded from canopy",
        "acronymn" => "unld",
        "unit" => "mm/tstep",
        "type" => MEAN
    ),
    "intcpt" => Dict(
        "longname" => "snow mass intercepted in the canopy",
        "acronymn" => "intc",
        "unit" => "mm/tstep",
        "type" => MEAN
    ),
    "Esrf" => Dict(
        "longname" => "moisture flux from the surface",
        "acronymn" => "esrf",
        "unit" => "kg/m^2/s",
        "type" => MEAN
    ),
    "Eveg" => Dict(
        "longname" => "moisture flux from vegetation",
        "acronymn" => "eveg",
        "unit" => "kg/m^2/s",
        "type" => MEAN
    ),
    "Gsoil" => Dict(
        "longname" => "heat flux into soil",
        "acronymn" => "ghsl",
        "unit" => "W/m^2",
        "type" => MEAN
    ),
    "Hsrf" => Dict(
        "longname" => "sensible heat flux from the surface",
        "acronymn" => "hesr",
        "unit" => "W/m^2",
        "type" => MEAN
    ),
    "LEsrf" => Dict(
        "longname" => "latent heat flux from the surface",
        "acronymn" => "lesr",
        "unit" => "W/m^2",
        "type" => MEAN
    ),
    "Rsrf" => Dict(
        "longname" => "net radiation at surface",
        "acronymn" => "rnsr",
        "unit" => "W/m^2",
        "type" => MEAN
    ),
    "Melt" => Dict(
        "longname" => "surface melt rate",
        "acronymn" => "emlt",
        "unit" => "kg/m^2/s",
        "type" => MEAN
    ),
    "SWsci" => Dict(
        "longname" => "subcanopy incoming shortwave radiation",
        "acronymn" => "swsc",
        "unit" => "W/m^2",
        "type" => MEAN
    ),
    "LWsci" => Dict(
        "longname" => "subcanopy incoming longwave radiation",
        "acronymn" => "lwsc",
        "unit" => "W/m^2",
        "type" => MEAN
    ),
    "LWveg" => Dict(
        "longname" => "net longwave radiation absorbed by vegetation",
        "acronymn" => "lwve",
        "unit" => "W/m^2",
        "type" => MEAN
    ),
    "Usc" => Dict(
        "longname" => "wind speed in canopy layer",
        "acronymn" => "uasc",
        "unit" => "m/s",
        "type" => MEAN
    ),
    "Sbveg" => Dict(
        "longname" => "sublimation from vegetation",
        "acronymn" => "sbve",
        "unit" => "mm/tstep",
        "type" => SUM
    ),
    "KH" => Dict(
        "longname" => "Eddy diffusivity for heat to the atmosphere",
        "acronymn" => "khag",
        "unit" => "m/s",
        "type" => MEAN
    ),
    "KHa" => Dict(
        "longname" => "Eddy diffusivity for heat from the canopy air space",
        "acronymn" => "khac",
        "unit" => "m/s",
        "type" => MEAN
    ),
    "KHg" => Dict(
        "longname" => "Eddy diffusivity for heat from the ground",
        "acronymn" => "khgr",
        "unit" => "m/s",
        "type" => MEAN
    ),
    "KHv" => Dict(
        "longname" => "Eddy diffusivity for heat from vegetation",
        "acronymn" => "khve",
        "unit" => "m/s",
        "type" => MEAN
    ),
    "KWg" => Dict(
        "longname" => "Eddy diffusivity for water from the ground",
        "acronymn" => "kwgr",
        "unit" => "m/s",
        "type" => MEAN
    ),
    "KWv" => Dict(
        "longname" => "Eddy diffusivity for water from vegetation",
        "acronymn" => "kwve",
        "unit" => "m/s",
        "type" => MEAN
    ),
    "snowdepthmin" => Dict(
        "longname" => "minimum snow depth at time step of swemin",
        "acronymn" => "hsmn",
        "unit" => "m",
        "type" => INSTANT
    ),
    "snowdepthmax" => Dict(
        "longname" => "maximum snow depth at time stemp of swemax",
        "acronymn" => "hsmx",
        "unit" => "m",
        "type" => INSTANT
    ),
    "snowdepthhist" => Dict(
        "longname" => "history of snow depth during last 14 days (most recent entries first)",
        "acronymn" => "hshs",
        "unit" => "m",
        "type" => INSTANT
    ),
    "swemin" => Dict(
        "longname" => "minimum swe during the season",
        "acronymn" => "swmn",
        "unit" => "mm",
        "type" => INSTANT
    ),
    "swemax" => Dict(
        "longname" => "maximum swe during the season",
        "acronymn" => "swmx",
        "unit" => "mm",
        "type" => INSTANT
    ),
    "swehist" => Dict(
        "longname" => "history of SWE during last 14 days",
        "acronymn" => "swhs",
        "unit" => "mm",
        "type" => INSTANT
    ),


    # HACK FOR TESTING SCF
    "fsnow_season_out" => Dict(
        "longname" => "fsnow_season_out",
        "acronymn" => "fsnow_season_out",
        "unit" => "fsnow_season_out",
        "type" => INSTANT
    ),
    "fsnow_nsnow_out" => Dict(
        "longname" => "fsnow_nsnow_out",
        "acronymn" => "fsnow_nsnow_out",
        "unit" => "fsnow_nsnow_out",
        "type" => INSTANT
    ),
    "fsnow_nsnow_recent_out" => Dict(
        "longname" => "fsnow_nsnow_recent_out",
        "acronymn" => "fsnow_nsnow_recent_out",
        "unit" => "fsnow_nsnow_recent_out",
        "type" => INSTANT
    )
    # HACK FOR TESTING SCF


)

# Variable access functions

get_variable_value(fsm::FSM, ::Val{:snowdepth}) = dropdims(sum(fsm.Ds, dims=1), dims=1) .* fsm.fsnow
get_variable_value(fsm::FSM, ::Val{:SWE}) = dropdims(sum(fsm.Sice .+ fsm.Sliq, dims=1), dims=1)
get_variable_value(fsm::FSM, ::Val{:Sliq_out}) = dropdims(sum(fsm.Sliq, dims=1), dims=1)
get_variable_value(fsm::FSM, ::Val{var}) where var = getfield(fsm, var)
get_variable_value(fsm::FSM, var::Symbol) = get_variable_value(fsm, Val(var))

get_variable_size(fsm::FSM, ::Val{:snowdepth}) = (fsm.Nx, fsm.Ny)
get_variable_size(fsm::FSM, ::Val{:SWE}) = (fsm.Nx, fsm.Ny)
get_variable_size(fsm::FSM, ::Val{:Sliq_out}) = (fsm.Nx, fsm.Ny)
get_variable_size(fsm::FSM, ::Val{var}) where var = size(getfield(fsm, var))
get_variable_size(fsm::FSM, var::Symbol) = get_variable_size(fsm, Val(var))

# Function for creating methods for accumulating and saving data

function make_saver(fsm::FSM, output_vars::Vector{String})

    invalid_vars = setdiff(output_vars, keys(AVAILABLE_OUTPUT_VARS))
    !isempty(invalid_vars) && error("Invalid variables: $(join(invalid_vars, ", "))")

    counter = 0

    enum_types = [AVAILABLE_OUTPUT_VARS[var]["type"] for var in output_vars]

    output_dicts = Dict[]
    for (var, type) in zip(output_vars, enum_types)
        size_var = get_variable_size(fsm, Symbol(var))
        output_dict = Dict(
            "data" => fill(0.0, size_var),
            "shortname" => var,
            "longname" => AVAILABLE_OUTPUT_VARS[var]["longname"],
            "acronymn" => AVAILABLE_OUTPUT_VARS[var]["acronymn"],
            "unit" => AVAILABLE_OUTPUT_VARS[var]["unit"],
            "type" => type == INSTANT ? "instant value at tend" :
                      type == SUM ? "sum from tstart to tend" :
                      "mean from tstart to tend"
        )
        push!(output_dicts, output_dict)
    end

    function accumulator(fsm::FSM)
        for (output_dict, enum_type) in zip(output_dicts, enum_types)
            value = get_variable_value(fsm, Symbol(output_dict["shortname"]))
            if enum_type == INSTANT
                output_dict["data"] .= value
            elseif enum_type == SUM || enum_type == MEAN
                output_dict["data"] .+= value
            end
        end
        counter += 1
        return nothing
    end

    function saver(fsm::FSM, t::DateTime, filepath::String; compress::Bool=true)
        file = nothing
        try
            file = matopen(filepath, "w"; compress=compress)
            for (output_dict, type) in zip(output_dicts, enum_types)
                type == MEAN && (output_dict["data"] .= output_dict["data"] ./ counter)
                var = output_dict["shortname"]
                write(file, var, output_dict)
            end
            write(file, "tend", Dates.format(t + Second(fsm.dt), "yyyy-mm-ddTHH:MM:SS"))
            write(file, "tstart", Dates.format(t - Second(fsm.dt * (counter - 1)), "yyyy-mm-ddTHH:MM:SS"))
        catch e
            @error "Failed to write output file $filepath: $e"
            rethrow(e)
        finally
            file !== nothing && close(file)
        end
        for output_dict in output_dicts
            output_dict["data"] .= 0.0
        end
        counter = 0
        return nothing
    end

    return accumulator, saver

end

function display_available_output_vars()

    descriptions = String[]
    variables = String[]
    units = String[]
    types = String[]

    for (var_name, var_info) in AVAILABLE_OUTPUT_VARS
        push!(descriptions, lowercase(var_info["longname"]))
        push!(variables, var_name)
        push!(units, var_info["unit"])
        push!(types, string(var_info["type"]))
    end

    perm = sortperm(descriptions)
    descriptions = descriptions[perm]
    variables = variables[perm]
    units = units[perm]
    types = types[perm]

    table = [descriptions variables units types]

    pretty_table(
        table,
        column_labels=["Description", "Variable", "Unit", "Type"],
        style=TextTableStyle(first_line_column_label=crayon"red bold");
        alignment=[:l, :l, :c, :c],
        display_size = (-1, -1)
    )

    println("\nTotal: $(length(variables)) variables available")

end
using ProgressMeter
using MAT
using Dates
using FSMOSHD

"""
    run_grid_simulation(; kwargs...)

Run a grid-based snow model simulation for open or forest tiles.

# Optional Arguments
- `settings`: Dict containing model configuration including "output_vars" for variables to save
  Default output_vars: ["snowdepth", "fsnow"]. See `AVAILABLE_OUTPUT_VARS` for all available options.
"""
function run_grid_simulation(;
    settings::Dict=Dict("tile" => "open", "config" => Dict("SNFRAC" => 0), "params" => Dict("wind_scaling" => 0.7), "output_vars" => ["snowdepth", "fsnow", "Roff", "Roff_snow"]),
    subfolder::String="default",
    base_folder::String="D:/julia",
    times::StepRange=DateTime(2024, 9, 1, 6):Hour(1):DateTime(2025, 6, 12, 6),
    Tf::Type=Float32,
    Ti::Type=Int32,
    verbose::Bool=true
)


    # Read landuse data
    landuse = prepare_landuse("K:/OSHD_AUX/DATA_LUS/OSHD_LUS_0250.mat")

    # Setup model
    Nx = size(landuse["dem"]["data"], 1)
    Ny = size(landuse["dem"]["data"], 2)

    fsm = setup(Tf, Ti, landuse, Nx, Ny, settings)

    met_curr = MET{Tf, Ti}(Nx=Nx, Ny=Ny)

    Sf24h = zeros(size(met_curr.Sf24h))
    Sf_history = zeros(size(met_curr.Sf_history))

    # Get output variables from settings
    output_vars = get(settings, "output_vars", ["snowdepth", "fsnow", "Roff", "Roff_snow"])
    
    # Create output directory
    mkpath(joinpath(base_folder, "FSM_HS_julia", subfolder))

    # Create accumulator and saver functions for storing model results
    accumulator, saver = make_saver(fsm, output_vars)

    # Run model
    for (i, t) in enumerate(times)

        elapsed_time = @elapsed begin

            # Prepare forcing data
            folder_icon = joinpath("K:/DATA_ICON/OUTPUT_OSHD_0250/PROCESSED_ANALYSIS/ICON_1EFA", Dates.format(t, "yyyy.mm"))
            filename_icon = searchdir(folder_icon, "ICONDATA_" * Dates.format(t, "yyyymmddHHMM"))
            met_single = matread(joinpath(folder_icon, filename_icon[1]))

            Sdir = met_single["sdri"]["data"]     # "direct shortwave radiation, per inclined surface area, within topography, above canopy"
            Sdif = met_single["sdfd"]["data"]     # "diffuse shortwave radiation, per horizontal surface area, within topography, above canopy"
            Sdird = met_single["sdrd"]["data"]    # "direct shortwave radiation, per horizontal surface area, within topography, above canopy"
            LW = met_single["lwrs"]["data"]       # "longwave radiation, above topography"
            Rf = met_single["prfc"]["data"]       # "rainfall"
            Sf = met_single["psfc"]["data"]       # "snowfall (snow + graupel)"
            Ta = met_single["tais"]["data"]       # "air temperature"
            RH = met_single["rhus"]["data"]       # "relative humidity"
            Ua = met_single["wnss"]["data"]       # "wind speed"
            Ps = met_single["pail"]["data"]       # "local air pressure"

            met_curr.Sdir[:, :] .= Sdir
            met_curr.Sdif[:, :] .= Sdif
            met_curr.Sdird[:, :] = Sdird
            met_curr.LW[:, :] .= LW
            met_curr.Sf[:, :] .= Sf
            met_curr.Rf[:, :] .= Rf
            met_curr.Ta[:, :] .= Ta
            met_curr.RH[:, :] .= RH
            met_curr.Ua[:, :] .= Ua
            met_curr.Ps[:, :] .= Ps

            if fsm.SNTRAN == 1
                if haskey(met_single,"wnsd")
                    Ua = met_single["wnsd"]["data"]       # "wind speed"
                else
                    Ua = met_single["wnss"]["data"]       # "wind speed"
                end
                if haskey(met_single,"wndd")
                    Udir = met_single["wndd"]["data"]     # "wind direction"
                else
                    Udir = met_single["wnds"]["data"]     # "wind direction"
                end
                met_curr.Ua[:, :] .= Ua
                met_curr.Udir[:, :] .= Udir
                Udir[isnan.(Udir)] .= 0
            end

            curr_hour = Dates.value(Hour(t)) + 1
            Sf24h .+= Sf
            Sf24h .-= Sf_history[:, :, curr_hour]
            Sf_history[:, :, curr_hour] = Sf
            met_curr.Sf24h[:, :] .= Sf24h

            # Forest-specific TVT data loading
            if settings["tile"] == "forest"
                folder_tvt = "K:/OSHD_AUX/DATA_CANRAD/OUTPUT_OSHD_0250/CR_2410_static/YYYY." * Dates.format(t, "mm")
                filename_tvt = searchdir(folder_tvt, "CANRAD_" * Dates.format(t, "mmddHHMM"))
                tvt_single = matread(joinpath(folder_tvt, filename_tvt[1]))
                met_curr.Tv[:, :] .= tvt_single["stdx"]["data"]
            end

            # Run model
            step!(fsm, met_curr, t)

            # Accumulate data for output
            accumulator(fsm)

            if hour(t) == 5
                output_filename = Dates.format(t + Dates.Hour(1), "yyyymmddHHMM") * "_output.mat"
                output_path = joinpath(base_folder, "FSM_HS_julia", subfolder, output_filename)
                saver(fsm, t, output_path)
            end

        end

        verbose && println(t, ", run time=", elapsed_time, " s")

    end

end
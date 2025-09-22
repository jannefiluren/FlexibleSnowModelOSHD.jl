using ProgressMeter
using MAT
using Dates
using FSMOSHD


const settings_default = Dict(
    "tile" => "open",
    "config" => Dict("SNFRAC" => 0),
    "params" => Dict("wind_scaling" => 0.7),
    "met_type" => :ICON,
    "met_folder" => "K:/DATA_ICON/OUTPUT_OSHD_0250/PROCESSED_ANALYSIS/ICON_1EFA",
    "met_prefix" => "ICON",
    "tvt_folder" => "K:/OSHD_AUX/DATA_CANRAD/OUTPUT_OSHD_0250/CR_2410_static",
    "lus_file" => "K:/OSHD_AUX/DATA_LUS/OSHD_LUS_0250.mat",
    "out_folder" => string(@__DIR__, "/../../FSM_HS_julia/default"),
    "output_vars" => ["snowdepth", "fsnow", "Roff", "meltflux_out"],
    )


"""
    run_grid_simulation(; kwargs...)

Run a grid-based snow model simulation for open, forest or glacier tiles.
"""
function run_grid_simulation(;
    settings::Dict=settings_default,
    times::StepRange=DateTime(2024, 9, 1, 6):Hour(1):DateTime(2025, 6, 12, 6),
    Tf::Type=Float32,
    Ti::Type=Int32,
    verbose::Bool=true,
)

    # Read landuse data
    landuse = prepare_landuse(settings["lus_file"])

    # Setup model
    Nx = size(landuse["dem"]["data"], 1)
    Ny = size(landuse["dem"]["data"], 2)

    fsm = setup(Tf, Ti, landuse, Nx, Ny, settings)
    met = MET{Tf,Ti}(Nx=Nx, Ny=Ny)

    # Get output variables from settings
    output_vars = get(settings, "output_vars", ["snowdepth", "fsnow", "Roff", "meltflux_out"])

    # Create output directory
    mkpath(settings["out_folder"])

    # Create accumulator and saver functions for storing model results
    accumulator, saver = make_saver(fsm, output_vars)

    # Run model
    for t in times

        elapsed_time = @elapsed begin

            # Read forcing data
            read_meteo!(met, fsm, t, settings)

            # Run model
            step!(fsm, met, t)

            # Accumulate data for output
            accumulator(fsm)

            # Save results
            if hour(t) == 5
                output_filename = Dates.format(t + Dates.Hour(1), "yyyymmddHHMM") * "_output.mat"
                output_filename = joinpath(settings["out_folder"], output_filename)
                saver(fsm, t, output_filename)
            end

        end

        verbose && println(t, ", run time=", elapsed_time, " s")

    end

end
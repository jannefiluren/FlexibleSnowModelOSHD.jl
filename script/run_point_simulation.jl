using FSMOSHD
using ProgressMeter
using MAT
using Dates
using CSV
using Tables

"""
    run_point_simulation(; kwargs...)

Run a point-based snow model simulation for station data with optional comparison to reference results.
"""
function run_point_simulation(;
    times::StepRange = DateTime(2024, 9, 1, 6):Hour(1):DateTime(2025, 6, 12, 6),
    landuse_file::String = "K:/OSHD_AUX/DATA_LUS/OSHD_LUS_STAT.mat",
    meteo_base_path::String = "K:/DATA_ICON/OUTPUT_OSHD_STAT/PROCESSED_ANALYSIS/ICON_1EFA",
    base_folder::String = "D:/julia",
    reference_path::String = "FSM_HS_matlab/LATEST_00h_RUN/OUTPUT_OSHD_STAT/RESULTS_01h_opn",
    write_outputs::Bool = true,
    Tf::Type=Float32,
    Ti::Type=Int32,
    verbose::Bool = true
)

    # Check if output path exist
    isdir(base_folder) || error("Base folder does not exist")

    # Prepare landuse data
    landuse = prepare_landuse(landuse_file)
    landuse_cropped = crop_landuse_to_domain(landuse)
    nstat = length(landuse_cropped["acro"])

    # Load meteorological data
    met_data = load_meteorological_data(times, meteo_base_path, nstat, verbose)

    # Setup model
    fsm = setup(Tf, Ti, landuse_cropped, nstat, 1, Dict("tile" => "open", "params" => Dict("wind_scaling" => 0.7)))
    met = MET{Tf, Ti}(Nx=nstat)

    # Run simulation
    snowdepth = run_simulation_loop(fsm, met, met_data, times, verbose)

    # Load reference results if requested
    snowdepth_ref = nothing
    if !isempty(reference_path)
        reference_path = joinpath(base_folder, reference_path)
        snowdepth_ref = load_reference_results(times, reference_path, nstat, verbose)
    end

    # Write outputs if requested
    if write_outputs
        write_simulation_outputs(snowdepth, snowdepth_ref, base_folder)
    end

    return (
        snowdepth = snowdepth,
        snowdepth_ref = snowdepth_ref,
        times = collect(times)
    )
end

"""
    load_meteorological_data(times, base_path, nstat, verbose=true)

Load meteorological forcing data for all time steps and stations.
"""
function load_meteorological_data(times, base_path, nstat, verbose=true)
    # Pre-allocate arrays
    Sdir = zeros(length(times), nstat)
    Sdif = zeros(length(times), nstat)
    Sdird = zeros(length(times), nstat)
    LW = zeros(length(times), nstat)
    Sf = zeros(length(times), nstat)
    Rf = zeros(length(times), nstat)
    Ta = zeros(length(times), nstat)
    RH = zeros(length(times), nstat)
    Ua = zeros(length(times), nstat)
    Ps = zeros(length(times), nstat)

    progress_bar = verbose ? Progress(length(times), desc="Loading meteorological data...") : nothing

    for (i, t) in enumerate(times)
        folder = joinpath(base_path, Dates.format(t, "yyyy.mm"))
        filename = searchdir(folder, "ICONDATA_" * Dates.format(t, "yyyymmddHHMM") * "_C1EFA_")
        
        met_single = matread(joinpath(folder, filename[1]))

        Sdir[i, :] = met_single["sdri"]["data"]
        Sdif[i, :] = met_single["sdfd"]["data"]
        Sdird[i, :] = met_single["sdrd"]["data"]
        LW[i, :] = met_single["lwrs"]["data"]
        Rf[i, :] = met_single["prfc"]["data"]
        Sf[i, :] = met_single["psfc"]["data"]
        Ta[i, :] = met_single["tais"]["data"]
        RH[i, :] = met_single["rhus"]["data"]
        Ua[i, :] = met_single["wnss"]["data"]
        Ps[i, :] = met_single["pail"]["data"]

        if verbose
            next!(progress_bar)
        end
    end

    return (
        Sdir = Sdir, Sdif = Sdif, Sdird = Sdird, LW = LW,
        Sf = Sf, Rf = Rf, Ta = Ta, RH = RH, Ua = Ua, Ps = Ps
    )
end

"""
    run_simulation_loop(fsm, met, met_data, times, verbose=true)

Execute the main simulation loop with snowfall tracking.
"""
function run_simulation_loop(fsm, met, met_data, times, verbose=true)
    # Setup snowfall tracking
    Sf24h = zeros(size(met.Sf24h))
    Sf_history = zeros(size(met.Sf_history))
    
    # Pre-allocate results
    snowdepth = zeros(length(times), size(fsm.Ds, 2))
    
    progress_bar = verbose ? Progress(length(times), desc="Running snow model...") : nothing

    for (i, t) in enumerate(times)
        # Set forcing data
        met.Sdir[:, :] .= met_data.Sdir[i, :]
        met.Sdif[:, :] .= met_data.Sdif[i, :]
        met.Sdird[:, :] = met_data.Sdird[i, :]
        met.LW[:, :] .= met_data.LW[i, :]
        met.Sf[:, :] .= met_data.Sf[i, :]
        met.Rf[:, :] .= met_data.Rf[i, :]
        met.Ta[:, :] .= met_data.Ta[i, :]
        met.RH[:, :] .= met_data.RH[i, :]
        met.Ua[:, :] .= met_data.Ua[i, :]
        met.Ps[:, :] .= met_data.Ps[i, :]

        # Update 24-hour snowfall tracking
        curr_hour = Dates.value(Hour(t)) + 1
        Sf24h .+= met_data.Sf[i, :]
        Sf24h .-= Sf_history[:, :, curr_hour]
        Sf_history[:, :, curr_hour] = met_data.Sf[i, :]
        met.Sf24h[:, :] .= Sf24h

        # Run model step
        step!(fsm, met, t)

        # Store results
        snowdepth[i, :] = dropdims(sum(fsm.Ds, dims=1), dims=3)

        if verbose
            next!(progress_bar)
        end
    end

    return snowdepth
end

"""
    load_reference_results(times, reference_path, nstat, verbose=true)

Load reference results for comparison.
"""
function load_reference_results(times, reference_path, nstat, verbose=true)
    snowdepth_ref = zeros(length(times), nstat)
    
    progress_bar = verbose ? Progress(length(times), desc="Loading reference results...") : nothing

    for (i, t) in enumerate(times)
        filename = searchdir(reference_path, "MODELDATA_" * Dates.format(t, "yyyymmddHHMM"))
        
        file = matopen(joinpath(reference_path, filename[1]))
        hsnt = read(file, "hsnt")
        close(file)

        snowdepth_ref[i, :] = hsnt["data"]

        if verbose
            next!(progress_bar)
        end
    end

    return snowdepth_ref
end

"""
    write_simulation_outputs(snowdepth, snowdepth_ref, base_folder)

Write simulation results to CSV files.
"""
function write_simulation_outputs(snowdepth, snowdepth_ref, base_folder)
    output_dir = joinpath(base_folder, "FSM_HS_matlab")
    mkpath(output_dir)

    CSV.write(joinpath(output_dir, "snowdepth_julia.csv"), Tables.table(snowdepth))
    
    if snowdepth_ref !== nothing
        CSV.write(joinpath(output_dir, "snowdepth_matlab.csv"), Tables.table(snowdepth_ref))
    end
end

"""
    compute_max_difference(snowdepth, snowdepth_ref; by_station=false)

Compute maximum difference between simulation and reference results.
"""
function compute_max_difference(snowdepth, snowdepth_ref; by_station::Bool=false)
    if snowdepth_ref === nothing
        error("Reference data is required to compute differences")
    end
    
    if size(snowdepth) != size(snowdepth_ref)
        error("Snowdepth arrays must have the same dimensions")
    end
    
    diff_array = abs.(snowdepth .- snowdepth_ref)
    
    if by_station
        # Return maximum difference for each station (across all times)
        return [maximum(diff_array[:, i]) for i in 1:size(diff_array, 2)]
    else
        # Return single maximum difference across all data
        return maximum(diff_array)
    end
end

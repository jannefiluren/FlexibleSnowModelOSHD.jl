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
  
Outputs the snowheight at 5 AM each day to .mat files in the specified output folder.
Optionally, a datetime object hs_out_time can be provided to return the snow height at that time and stop the simulation.

"""
function run_grid_simulation(;
    settings::Dict=Dict("tile" => "open"),
    out_folder::String="D:/snowinflow_project/snowinflow_data",
    times::StepRange=DateTime(2024, 9, 1, 6):Hour(1):DateTime(2025, 6, 12, 6),
    Tf::Type=Float32,
    Ti::Type=Int32,
    verbose::Bool=true,
    daily_output::Bool=true,
    lus_file::String,
    met_folder::String,
    hs_out_times::Vector{DateTime}=[],
)


    # Read landuse data
    landuse = prepare_landuse(lus_file)

    # Setup model
    Nx = size(landuse["dem"]["data"], 1)
    Ny = size(landuse["dem"]["data"], 2)

    fsm = setup(Tf, Ti, landuse, Nx, Ny, settings)

    met_curr = MET{Tf, Ti}(Nx=Nx, Ny=Ny)

    Sf24h = zeros(size(met_curr.Sf24h))
    Sf_history = zeros(size(met_curr.Sf_history))

    # Create output directory
    mkpath(out_folder)
     output_vars = get(settings, "output_vars", ["snowdepth", "fsnow", "Roff", "meltflux_out"])
    accumulator, saver = make_saver(fsm, output_vars)

    # Preallocate output array
    hs_out = Array{Tf}(undef, Nx, Ny, length(hs_out_times))

    #timers
    tot_time_model = 0.0
    tot_time_met = 0.0
    tot_time_output = 0.0
    tot_time = 0.0

    # Run model
    println("Starting model run...")
    for (i, t) in enumerate(times)

        elapsed_time_met = @elapsed begin

            # Prepare forcing data
            folder_icon = joinpath(met_folder, "PROCESSED_ANALYSIS", Dates.format(t, "yyyy.mm"))
            filename_icon = searchdir(folder_icon, "METNORDICDATA_" * Dates.format(t, "yyyymmddHHMM"))
            met_single = matread(joinpath(folder_icon, filename_icon[1]))

            Sdir = round.(met_single["sdrd"]["data"], digits=2)     # "direct shortwave radiation, per inclined surface area, within topography, above canopy"
            Sdif = round.(met_single["sdfd"]["data"], digits=2)     # "diffuse shortwave radiation, per horizontal surface area, within topography, above canopy"
            LW = round.(met_single["lwrs"]["data"], digits=2)       # "longwave radiation, above topography"
            Ptot = round.(met_single["prcs"]["data"], digits=2)     # "precipitation"
            Ta = round.(met_single["tais"]["data"], digits=2)       # "air temperature"
            RH = round.(met_single["rhus"]["data"], digits=2)       # "relative humidity"
            Ua = round.(met_single["wnss"]["data"], digits=2)       # "wind speed"
            Ps = round.(met_single["pail"]["data"], digits=2)       # "local air pressure"

            thres_prec = 1.04 + 273.15
            m_prec = 0.15
            Sf = 1.0 ./ (1 .+ exp.((Ta .- thres_prec) ./ m_prec)) .* Ptot;
            Rf = Ptot .- Sf;
            
            Sf .*= settings["sfscaling"] #apply scaling factor
            # print(settings["sfscaling"])
            
            met_curr.Sdir[:, :] .= Sdir
            met_curr.Sdif[:, :] .= Sdif
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
        end

        elapsed_time_model = @elapsed begin
            # Run model
            step!(fsm, met_curr, t)
        end
        
        elapsed_time_output = @elapsed begin
            # Output data at 5 AM
            if hour(t) == 5 && daily_output
                output_filename = Dates.format(t + Dates.Hour(1), "yyyymmddHHMM") * "_output.mat"
                output_path = joinpath(out_folder, output_filename)
                saver(fsm, t, output_path)
            end
            if t in hs_out_times
                #save in the output array
                hs_out[:,:,findfirst(t .== hs_out_times)] = dropdims(sum(fsm.Ds, dims=1), dims=1) .* fsm.fsnow
            end
        end
        loop_time = elapsed_time_met+elapsed_time_model+elapsed_time_output
        verbose && print("\r\e[2K") #clear line
        verbose && print(t, ", run time=", round(loop_time; digits=3), " s")#, "(meteo: ", round(elapsed_time_met/tot_time*100), "%, model: ", round(elapsed_time_model/tot_time*100), "%, output: ", round(elapsed_time_output/tot_time*100), "%)")
        tot_time_met += elapsed_time_met
        tot_time_model += elapsed_time_model
        tot_time_output += elapsed_time_output
        tot_time += loop_time
    end
    verbose && print("\r\e[2K") #clear line
    verbose && println("Finished! Total run time: ", round(tot_time; digits=3), " s (meteo: ", round(tot_time_met/tot_time*100), "%, model: ", round(tot_time_model/tot_time*100), "%, output: ", round(tot_time_output/tot_time*100), "%)")

    return hs_out
end
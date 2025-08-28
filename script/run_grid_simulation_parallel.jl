using MAT
using Dates
using FSMOSHD
using ChunkSplitters

"""
    create_fsm_chunk(fsm_full::FSM{Tf, Ti}, chunk_range::UnitRange{Int}) where {Tf, Ti}

Create FSM chunk for spatial domain decomposition.
"""
function create_fsm_chunk(fsm_full::FSM{Tf, Ti}, chunk_range::UnitRange{Int}) where {Tf, Ti}
    fsm_chunk = deepcopy(fsm_full)
    fsm_chunk.Nx = length(chunk_range)
    fsm_chunk.Ny = fsm_full.Ny
    
    for field in fieldnames(typeof(fsm_full))
        val = getfield(fsm_full, field)
        if val isa Array
            if ndims(val) == 2 && size(val, 1) == fsm_full.Nx && size(val, 2) == fsm_full.Ny
                setfield!(fsm_chunk, field, copy(val[chunk_range, :]))
            elseif ndims(val) == 3 && size(val, 2) == fsm_full.Nx && size(val, 3) == fsm_full.Ny
                setfield!(fsm_chunk, field, copy(val[:, chunk_range, :]))
            end
        end
    end
    
    return fsm_chunk
end

"""
    create_met_chunk(chunk_Nx, Ny, ::Type{Tf}, ::Type{Ti}) where {Tf, Ti}

Create empty MET structure for chunk.
"""
function create_met_chunk(chunk_Nx, Ny, ::Type{Tf}, ::Type{Ti}) where {Tf, Ti}
    return MET{Tf, Ti}(Nx=chunk_Nx, Ny=Ny)
end

"""
    fill_met_chunk!(met_chunk::MET{Tf, Ti}, met_full::MET{Tf, Ti}, chunk_range::UnitRange{Int}) where {Tf, Ti}

Copy meteorological data to chunk.
"""
function fill_met_chunk!(met_chunk::MET{Tf, Ti}, met_full::MET{Tf, Ti}, chunk_range::UnitRange{Int}) where {Tf, Ti}
    @views begin
        met_chunk.LW .= met_full.LW[chunk_range, :]
        met_chunk.Ps .= met_full.Ps[chunk_range, :]
        met_chunk.Qa .= met_full.Qa[chunk_range, :]
        met_chunk.Rf .= met_full.Rf[chunk_range, :]
        met_chunk.RH .= met_full.RH[chunk_range, :]
        met_chunk.Sdif .= met_full.Sdif[chunk_range, :]
        met_chunk.Sdir .= met_full.Sdir[chunk_range, :]
        met_chunk.Sdird .= met_full.Sdird[chunk_range, :]
        met_chunk.Sf .= met_full.Sf[chunk_range, :]
        met_chunk.Sf24h .= met_full.Sf24h[chunk_range, :]
        met_chunk.Ta .= met_full.Ta[chunk_range, :]
        met_chunk.Tv .= met_full.Tv[chunk_range, :]
        met_chunk.Ua .= met_full.Ua[chunk_range, :]
        met_chunk.Sf_history .= met_full.Sf_history[chunk_range, :, :]
    end
end

"""
    extract_outputs!(hs_output::Array{Tf,2}, fsnow_output::Array{Tf,2}, 
                    fsm_chunks::Vector{FSM{Tf, Ti}}, chunk_ranges) where {Tf, Ti}

Extract simulation outputs from chunks.
"""
function extract_outputs!(hs_output::Array{Tf,2}, fsnow_output::Array{Tf,2}, 
                         fsm_chunks::Vector{FSM{Tf, Ti}}, chunk_ranges) where {Tf, Ti}
    for (chunk_idx, chunk_range) in enumerate(chunk_ranges)
        fsm_chunk = fsm_chunks[chunk_idx]
        
        fsnow_output[chunk_range, :] .= fsm_chunk.fsnow
        @views hs_output[chunk_range, :] .= dropdims(sum(fsm_chunk.Ds, dims=1), dims=1) .* fsm_chunk.fsnow
    end
end

"""
    run_grid_simulation_parallel(; kwargs...)

Run a grid-based snow model simulation with parallel execution using spatial domain decomposition.
"""
function run_grid_simulation_parallel(;
    times::StepRange=DateTime(2024, 9, 1, 6):Hour(1):DateTime(2025, 6, 12, 6),
    settings::Dict=Dict("tile" => "open", "config" => Dict("SNFRAC" => 0), "params" => Dict("wind_scaling" => 0.7)),
    landuse_file::String="K:/OSHD_AUX/DATA_LUS/OSHD_LUS_0250.mat",
    base_folder::String="D:/julia",
    subfolder::String="parallel_$(Dates.format(now(), "yyyymmdd_HHMMSS"))",
    Tf::Type=Float32,
    Ti::Type=Int32,
    verbose::Bool=true,
    nchunks::Int=Threads.nthreads()
)


    if verbose
        println("Starting parallel simulation with $(Threads.nthreads()) threads")
    end

    # Read landuse data
    landuse = prepare_landuse(landuse_file)

    # Setup model
    Nx = size(landuse["dem"]["data"], 1)
    Ny = size(landuse["dem"]["data"], 2)

    if verbose
        println("Domain size: $Nx × $Ny grid cells")
    end

    fsm = setup(Tf, Ti, landuse, Nx, Ny, settings)
    met_curr = MET{Tf, Ti}(Nx=Nx, Ny=Ny)

    # Initialize snowfall tracking arrays
    Sf24h = zeros(size(met_curr.Sf24h))
    Sf_history = zeros(size(met_curr.Sf_history))

    # Calculate domain partitioning
    chunk_ranges = index_chunks(1:Nx; n=nchunks)
    
    if verbose
        println("Domain partitioned into $(length(chunk_ranges)) chunks:")
        for (i, range) in enumerate(chunk_ranges)
            println("  Chunk $i: indices $(first(range))-$(last(range)) ($(length(range)) cells)")
        end
    end

    # Create structs for domain chunks
    if verbose
        println("Creating $(length(chunk_ranges)) FSM and MET chunks...")
    end
    
    fsm_chunks = Vector{FSM{Tf, Ti}}(undef, length(chunk_ranges))
    met_chunks = Vector{MET{Tf, Ti}}(undef, length(chunk_ranges))
    
    for (chunk_idx, chunk_range) in enumerate(chunk_ranges)
        chunk_Nx = length(chunk_range)
        
        fsm_chunks[chunk_idx] = create_fsm_chunk(fsm, chunk_range)
        
        met_chunks[chunk_idx] = create_met_chunk(chunk_Nx, Ny, Tf, Ti)
        
        if verbose
            println("  ✓ Chunk $chunk_idx: $(chunk_Nx) × $Ny cells")
        end
    end
    
    # Pre-allocate output arrays
    hs_output = zeros(Tf, Nx, Ny)
    fsnow_output = zeros(Tf, Nx, Ny)
    
    # Constants for meteorological data loading
    icon_base = "K:/DATA_ICON/OUTPUT_OSHD_0250/PROCESSED_ANALYSIS/ICON_1EFA"
    tvt_base = "K:/OSHD_AUX/DATA_CANRAD/OUTPUT_OSHD_0250/CR_2410_static/YYYY."
    output_path = joinpath(base_folder, "FSM_HS_julia", subfolder)
    
    met_var_mapping = [
            ("sdri", met_curr.Sdir), ("sdfd", met_curr.Sdif), ("sdrd", met_curr.Sdird),
            ("lwrs", met_curr.LW), ("prfc", met_curr.Rf), ("psfc", met_curr.Sf),
            ("tais", met_curr.Ta), ("rhus", met_curr.RH), ("wnss", met_curr.Ua), ("pail", met_curr.Ps)
        ]

    # Create output directory
    mkpath(output_path)
        
    # Time loop
    if verbose
        println("Starting time loop with parallel execution...")
    end
    
    for (i, t) in enumerate(times)

        elapsed_total = @elapsed begin

            # Load forcing data
            elapsed_loading = @elapsed begin
                
                folder_icon = joinpath(icon_base, Dates.format(t, "yyyy.mm"))
                filename_icon = searchdir(folder_icon, "ICONDATA_" * Dates.format(t, "yyyymmddHHMM"))
                mat_filepath = joinpath(folder_icon, filename_icon[1])
                
                Threads.@threads for met_var in met_var_mapping
                    file = matopen(mat_filepath)
                    var_name, target_array = met_var
                    data = read(file, var_name)["data"]
                    close(file)
                    target_array[:, :] .= data
                end
    
                curr_hour = Dates.value(Hour(t)) + 1
                Sf24h .+= met_curr.Sf
                Sf24h .-= Sf_history[:, :, curr_hour]
                Sf_history[:, :, curr_hour] = met_curr.Sf
                met_curr.Sf24h[:, :] .= Sf24h
    
                if settings["tile"] == "forest"
                    folder_tvt = tvt_base * Dates.format(t, "mm")
                    filename_tvt = searchdir(folder_tvt, "CANRAD_" * Dates.format(t, "mmddHHMM"))
                    tvt_single = matread(joinpath(folder_tvt, filename_tvt[1]))
                    met_curr.Tv[:, :] .= tvt_single["stdx"]["data"]
                end
    
                for (chunk_idx, chunk_range) in enumerate(chunk_ranges)
                    fill_met_chunk!(met_chunks[chunk_idx], met_curr, chunk_range)
                end

            end

            # Run model
            elapsed_simulation = @elapsed begin
                
                Threads.@threads for chunk_idx in 1:length(chunk_ranges)
                    step!(fsm_chunks[chunk_idx], met_chunks[chunk_idx], t)
                end

            end

            # Save outputs
            elapsed_saving = @elapsed begin
                
                if hour(t) == 5
                    extract_outputs!(hs_output, fsnow_output, fsm_chunks, chunk_ranges)
                    filename = Dates.format(t + Hour(1), "yyyymmddHHMM") * "_output.mat"
                    matwrite(joinpath(output_path, filename), 
                        Dict("hs" => hs_output, "fsnow" => fsnow_output); compress=true)
                end

            end

        end

        if verbose
            println("$t, run time=$(round(elapsed_total, digits=2)) s ($(length(chunk_ranges)) chunks) ✓")
            println("Components, loading=$(round(elapsed_loading, digits=2)), simulation=$(round(elapsed_simulation, digits=2)), saving=$(round(elapsed_saving, digits=2))")
        end

    end

    if verbose
        println("Parallel simulation completed. Output saved to: $output_path")
    end

end
using MAT
using Dates
using FSMOSHD

"""
    partition_domain(Nx::Int, nchunks::Int) -> Vector{UnitRange{Int}}

Partition a domain of size Nx into nchunks roughly equal chunks.
"""
function partition_domain(Nx::Int, nchunks::Int)
    if nchunks >= Nx
        return [i:i for i in 1:Nx]
    end
    
    chunk_size = Nx ÷ nchunks
    remainder = Nx % nchunks
    
    ranges = Vector{UnitRange{Int}}()
    start_idx = 1
    
    for i in 1:nchunks
        extra = i <= remainder ? 1 : 0
        end_idx = start_idx + chunk_size - 1 + extra
        push!(ranges, start_idx:end_idx)
        start_idx = end_idx + 1
    end
    
    return ranges
end

"""
    create_fsm_chunk_safe(fsm_full::FSM{Tf,Ti}, i_range::UnitRange{Int}) where {Tf,Ti}

Create an FSM chunk with ALL required fields properly initialized.
Uses a safer approach that copies EVERYTHING from the original FSM.
"""
function create_fsm_chunk_safe(fsm_full::FSM{Tf,Ti}, i_range::UnitRange{Int}) where {Tf,Ti}
    chunk_Nx = length(i_range)
    Ny = fsm_full.Ny
    
    # Start with a complete copy of the FSM and then adjust dimensions and arrays
    fsm_chunk = deepcopy(fsm_full)
    
    # Update dimensions
    fsm_chunk.Nx = chunk_Nx
    fsm_chunk.Ny = Ny
    
    # Replace all 2D arrays (Nx, Ny) with chunks  
    fsm_chunk.fsnow = copy(fsm_full.fsnow[i_range, :])
    fsm_chunk.alb = copy(fsm_full.alb[i_range, :])
    fsm_chunk.Tsrf = copy(fsm_full.Tsrf[i_range, :])
    fsm_chunk.Nsnow = copy(fsm_full.Nsnow[i_range, :])
    
    # Replace all 3D arrays (Nsmax, Nx, Ny) with chunks
    fsm_chunk.Ds = copy(fsm_full.Ds[:, i_range, :])
    fsm_chunk.Sice = copy(fsm_full.Sice[:, i_range, :])
    fsm_chunk.Sliq = copy(fsm_full.Sliq[:, i_range, :])
    fsm_chunk.Tsnow = copy(fsm_full.Tsnow[:, i_range, :])
    fsm_chunk.rgrn = copy(fsm_full.rgrn[:, i_range, :])
    fsm_chunk.histowet = copy(fsm_full.histowet[:, i_range, :])
    fsm_chunk.ksnow = copy(fsm_full.ksnow[:, i_range, :])
    
    # Replace all 3D arrays (Nsoil, Nx, Ny) with chunks
    fsm_chunk.Tsoil = copy(fsm_full.Tsoil[:, i_range, :])
    fsm_chunk.theta = copy(fsm_full.theta[:, i_range, :])
    fsm_chunk.csoil = copy(fsm_full.csoil[:, i_range, :])
    fsm_chunk.ksoil = copy(fsm_full.ksoil[:, i_range, :])
    
    # Replace all other 2D arrays with chunks
    for field in fieldnames(typeof(fsm_full))
        val = getfield(fsm_full, field)
        if val isa Array && ndims(val) == 2 && size(val, 1) == fsm_full.Nx && size(val, 2) == fsm_full.Ny
            # This is a 2D array with dimensions (Nx, Ny) - create chunk
            setfield!(fsm_chunk, field, copy(val[i_range, :]))
        end
    end
    
    return fsm_chunk
end

"""
    create_met_chunk(chunk_Nx::Ti, Ny::Ti, ::Type{Tf}, ::Type{Ti}) where {Tf,Ti}

Create an empty MET structure for a chunk with the specified dimensions.
"""
function create_met_chunk(chunk_Nx::Ti, Ny::Ti, ::Type{Tf}, ::Type{Ti}) where {Tf,Ti}
    return MET{Tf,Ti}(Nx=chunk_Nx, Ny=Ny)
end

"""
    fill_met_chunk!(met_chunk::MET{Tf,Ti}, met_full::MET{Tf,Ti}, i_range::UnitRange{Int}) where {Tf,Ti}

Copy meteorological data from the full MET structure to a chunk.
"""
function fill_met_chunk!(met_chunk::MET{Tf,Ti}, met_full::MET{Tf,Ti}, i_range::UnitRange{Int}) where {Tf,Ti}
    # Copy 2D meteorological fields
    met_chunk.LW .= met_full.LW[i_range, :]
    met_chunk.Ps .= met_full.Ps[i_range, :]
    met_chunk.Qa .= met_full.Qa[i_range, :]
    met_chunk.Rf .= met_full.Rf[i_range, :]
    met_chunk.RH .= met_full.RH[i_range, :]
    met_chunk.Sdif .= met_full.Sdif[i_range, :]
    met_chunk.Sdir .= met_full.Sdir[i_range, :]
    met_chunk.Sdird .= met_full.Sdird[i_range, :]
    met_chunk.Sf .= met_full.Sf[i_range, :]
    met_chunk.Sf24h .= met_full.Sf24h[i_range, :]
    met_chunk.Ta .= met_full.Ta[i_range, :]
    met_chunk.Tv .= met_full.Tv[i_range, :]
    met_chunk.Ua .= met_full.Ua[i_range, :]
    
    # Copy 3D field (Sf_history)
    met_chunk.Sf_history .= met_full.Sf_history[i_range, :, :]
    
    return nothing
end

"""
    copy_results_back_safe!(fsm_full::FSM{Tf,Ti}, fsm_chunk::FSM{Tf,Ti}, i_range::UnitRange{Int}) where {Tf,Ti}

Safely copy results from a chunk back to the full FSM structure.
"""
function copy_results_back_safe!(fsm_full::FSM{Tf,Ti}, fsm_chunk::FSM{Tf,Ti}, i_range::UnitRange{Int}) where {Tf,Ti}
    # Copy back all 2D arrays that might have changed
    for field in fieldnames(typeof(fsm_full))
        val_full = getfield(fsm_full, field)
        if val_full isa Array && ndims(val_full) == 2 && size(val_full, 1) == fsm_full.Nx && size(val_full, 2) == fsm_full.Ny
            val_chunk = getfield(fsm_chunk, field)
            val_full[i_range, :] .= val_chunk
        end
    end
    
    # Copy back 3D snow arrays (Nsmax, Nx, Ny)
    fsm_full.Ds[:, i_range, :] .= fsm_chunk.Ds
    fsm_full.Sice[:, i_range, :] .= fsm_chunk.Sice
    fsm_full.Sliq[:, i_range, :] .= fsm_chunk.Sliq
    fsm_full.Tsnow[:, i_range, :] .= fsm_chunk.Tsnow
    fsm_full.rgrn[:, i_range, :] .= fsm_chunk.rgrn
    fsm_full.histowet[:, i_range, :] .= fsm_chunk.histowet
    fsm_full.ksnow[:, i_range, :] .= fsm_chunk.ksnow
    
    # Copy back 3D soil arrays (Nsoil, Nx, Ny)
    fsm_full.Tsoil[:, i_range, :] .= fsm_chunk.Tsoil
    fsm_full.theta[:, i_range, :] .= fsm_chunk.theta
    fsm_full.csoil[:, i_range, :] .= fsm_chunk.csoil
    fsm_full.ksoil[:, i_range, :] .= fsm_chunk.ksoil
    
    return nothing
end

"""
    run_grid_simulation_parallel(; kwargs...)

Run a grid-based snow model simulation with parallel execution using spatial domain decomposition.

This function splits the spatial domain into chunks and processes each chunk in parallel
using different threads. Each thread operates on a contiguous block of grid cells,
providing significant speedup for large domains.

# Arguments
- `settings::Dict`: Configuration dictionary
- `subfolder::String`: Output subfolder (auto-generated if empty)  
- `base_folder::String="D:/julia"`: Base output folder
- `times::StepRange`: Time range for simulation
- `Tf::Type=Float32`: Numerical precision for floats
- `Ti::Type=Int32`: Numerical precision for integers
- `verbose::Bool=true`: Display progress information
- `nchunks::Int=Threads.nthreads()`: Number of spatial chunks to use

# Performance Notes
- Best performance typically achieved when `nchunks` equals the number of CPU cores
- Memory usage scales with number of chunks (each chunk gets full FSM copy)
- Speedup depends on domain size - larger domains benefit more from parallelization

# Example
```julia
# Run with 4 threads and 4 chunks
julia -t 4 -e '
include("script/run_grid_simulation_parallel.jl");
run_grid_simulation_parallel(
    settings=Dict("tile" => "open"),
    times=DateTime(2024,12,1,6):Hour(1):DateTime(2024,12,2,6),
    nchunks=4
)'
```
"""
function run_grid_simulation_parallel(;
    settings::Dict=Dict("tile" => "open", "config" => Dict("SNFRAC" => 0)),
    subfolder::String="parallel_$(Dates.format(now(), "yyyymmdd_HHMMSS"))",
    base_folder::String="D:/julia",
    times::StepRange=DateTime(2024, 9, 1, 6):Hour(1):DateTime(2025, 6, 12, 6),
    Tf::Type=Float32,
    Ti::Type=Int32,
    verbose::Bool=true,
    nchunks::Int=Threads.nthreads()  # Default to match thread count
)

    # Helper function
    searchdir(path, key) = filter(x -> occursin(key, x), readdir(path))

    if verbose
        println("Starting parallel simulation with $(Threads.nthreads()) threads and $nchunks chunks")
    end

    # Read landuse data
    landuse = prepare_landuse("K:/OSHD_AUX/DATA_LUS/OSHD_LUS_0250.mat")

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
    chunk_ranges = partition_domain(Nx, nchunks)
    
    if verbose
        println("Domain partitioned into $(length(chunk_ranges)) chunks:")
        for (i, range) in enumerate(chunk_ranges)
            println("  Chunk $i: indices $(first(range))-$(last(range)) ($(length(range)) cells)")
        end
    end

    # PRE-CREATE CHUNKS using safe method
    if verbose
        println("Creating $(length(chunk_ranges)) FSM and MET chunks...")
    end
    
    fsm_chunks = Vector{FSM{Tf,Ti}}(undef, length(chunk_ranges))
    met_chunks = Vector{MET{Tf,Ti}}(undef, length(chunk_ranges))
    
    for (chunk_idx, i_range) in enumerate(chunk_ranges)
        chunk_Nx = length(i_range)
        
        # Create FSM chunk with SAFE method (deepcopy approach)
        fsm_chunks[chunk_idx] = create_fsm_chunk_safe(fsm, i_range)
        
        # Create empty MET chunk
        met_chunks[chunk_idx] = create_met_chunk(Ti(chunk_Nx), Ti(Ny), Tf, Ti)
        
        if verbose
            println("  ✓ Chunk $chunk_idx: $(chunk_Nx) × $Ny cells")
        end
    end

    # Create output directory
    mkpath(joinpath(base_folder, "FSM_HS_julia", subfolder))

    # TIME LOOP with enhanced error handling
    if verbose
        println("Starting time loop with parallel execution...")
    end
    
    for (i, t) in enumerate(times)

        elapsed_time = @elapsed begin

            # Load meteorological data
            folder_icon = joinpath("K:/DATA_ICON/OUTPUT_OSHD_0250/PROCESSED_ANALYSIS/ICON_1EFA", Dates.format(t, "yyyy.mm"))
            filename_icon = searchdir(folder_icon, "ICONDATA_" * Dates.format(t, "yyyymmddHHMM") * "_C1EFA_")
            met_single = matread(joinpath(folder_icon, filename_icon[1]))

            # ... (same meteorological data loading as before)
            Sdir = met_single["sdri"]["data"]
            Sdif = met_single["sdfd"]["data"]
            Sdird = met_single["sdrd"]["data"]
            LW = met_single["lwrs"]["data"]
            Rf = met_single["prfc"]["data"]
            Sf = met_single["psfc"]["data"]
            Ta = met_single["tais"]["data"]
            RH = met_single["rhus"]["data"]
            Ua = met_single["wnss"]["data"]
            Ps = met_single["pail"]["data"]

            # Update full meteorological arrays
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

            # Update 24-hour snowfall tracking
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

            # FILL METEOROLOGICAL DATA INTO CHUNKS
            for (chunk_idx, i_range) in enumerate(chunk_ranges)
                fill_met_chunk!(met_chunks[chunk_idx], met_curr, i_range)
            end

            # PARALLEL EXECUTION with try-catch for better error reporting
            elapsed_step = @elapsed begin
                Threads.@threads for chunk_idx in 1:length(chunk_ranges)
                    step!(fsm_chunks[chunk_idx], met_chunks[chunk_idx], t)
                end
            end
            println(elapsed_step)

            # COPY RESULTS BACK to main FSM
            for (chunk_idx, i_range) in enumerate(chunk_ranges)
                copy_results_back_safe!(fsm, fsm_chunks[chunk_idx], i_range)
            end

            # Output data at 5 AM
            if hour(t) == 5
                matwrite(joinpath(base_folder, "FSM_HS_julia", subfolder, Dates.format(t + Dates.Hour(1), "yyyymmddHHMM") * "_output.mat"),
                    Dict(
                        "hs" => dropdims(sum(fsm.Ds, dims=1), dims=1) .* fsm.fsnow,
                        "fsnow" => fsm.fsnow
                    ); compress=true)
            end

        end

        if verbose
            println("$t, run time=$(round(elapsed_time, digits=3)) s ($(length(chunk_ranges)) chunks) ✓")
        end

    end

    if verbose
        println("Parallel simulation completed. Output saved to: $(joinpath(base_folder, "FSM_HS_julia", subfolder))")
    end

end
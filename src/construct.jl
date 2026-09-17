# Fetch a required landuse field, with a clear error naming the missing key.
function require_field(landuse::Dict, key)
    haskey(landuse, key) || error("landuse is missing required field \"$key\"")
    return landuse[key]["data"]
end

"""
    Grid(Tf, landuse::Dict; kwargs...)

Build a `Grid` for a landuse domain; `Nx`, `Ny` are taken from the elevation array.
"""
function Grid(::Type{Tf}, landuse::Dict; kwargs...) where {Tf}
    Nx, Ny = size(require_field(landuse, "elevation"))
    return Grid{Tf, Vector{Tf}}(; kwargs..., Nx, Ny)
end

# Build the Surface struct from the landuse data. Which landuse fraction to load and whether to read
# canopy — is inferred from the scheme types.
function build_surface(grid::Grid{Tf}, landuse::Dict, land_cover, substrate, params) where {Tf}
    c = get_constants(Tf)
    sf = Surface{typeof(grid), Matrix{Tf}, Matrix{Float64}}(; grid = grid)

    # Terrain properties
    sf.fsky_terr .= Tf.(require_field(landuse, "skyvf"))
    sf.dem .= Tf.(require_field(landuse, "elevation"))
    sf.prec_multi .= require_field(landuse, "prec_multi")
    sf.slopemu .= Tf.(require_field(landuse, "slopemu"))
    sf.xi .= Tf.(require_field(landuse, "xi"))
    sf.Ld .= Tf.(require_field(landuse, "Ld"))

    # Tile fraction for the chosen tile (open keeps the default of ones)
    if land_cover isa ForestCover
        sf.tilefrac .= Tf.(require_field(landuse, "forest"))
    elseif substrate isa IceSubstrate
        sf.tilefrac .= Tf.(require_field(landuse, "glacier"))
    end

    # Canopy inputs (forest tile only)
    if land_cover isa ForestCover
        sf.fveg .= Tf.(require_field(landuse, "fveg"))
        sf.hcan .= Tf.(require_field(landuse, "hcan"))
        sf.lai .= Tf.(require_field(landuse, "lai"))
        sf.vfhp .= Tf.(require_field(landuse, "vfhp"))
        sf.fves .= Tf.(require_field(landuse, "fves"))
    end

    # Derived soil parameters
    mask = sf.fcly .+ sf.fsnd .> Tf(1)
    sf.fcly[mask] .= Tf(1) .- sf.fsnd[mask]
    sf.b .= Tf(3.1) .+ Tf(15.7) .* sf.fcly .- Tf(0.3) .* sf.fsnd
    sf.hcap_soil .= (Tf(2.128) .* sf.fcly .+ Tf(2.385) .* sf.fsnd) .* Tf(1.0e6) ./ (sf.fcly .+ sf.fsnd)
    sf.sathh .= Tf(10) .^ (Tf(0.17) .- Tf(0.63) .* sf.fcly .- Tf(1.58) .* sf.fsnd)
    sf.Vsat .= Tf(0.505) .- Tf(0.037) .* sf.fcly .- Tf(0.142) .* sf.fsnd
    sf.Vcrit .= sf.Vsat .* (sf.sathh ./ Tf(3.364)) .^ (Tf(1) ./ sf.b)
    hcon_min = (c.hcon_clay .^ sf.fcly) .* (c.hcon_sand .^ (Tf(1) .- sf.fcly))
    sf.hcon_soil .= (c.hcon_air .^ sf.Vsat) .* (hcon_min .^ (Tf(1) .- sf.Vsat))

    # Derived canopy fields (forest tile only)
    if land_cover isa ForestCover
        forest = require_field(landuse, "forest")
        prec_multi = require_field(landuse, "prec_multi")
        sf.pmultf .= Tf.((1 .- (1 .- prec_multi) .* (1 .- forest * params.pmultf_for)) ./ prec_multi)
        sf.VAI[:, :] = sf.lai[:, :]
        sf.trcn[:, :] = Tf(1) .- Tf(0.9) .* sf.fveg[:, :]
        sf.fsky .= sf.vfhp ./ sf.trcn
        cmask = sf.fsky .> Tf(1)
        sf.trcn[cmask] .= sf.vfhp[cmask]
        sf.fsky[cmask] .= Tf(1)
        sf.canh[:, :] = Tf(12500) .* sf.VAI[:, :]
        sf.scap[:, :] = params.cvai .* sf.VAI[:, :]
    end

    return sf
end

# Validate the built Surface sturct against the chosen physics
function validate_surface(surface, params, land_cover)
    Tf = eltype(surface.dem)

    if land_cover isa ForestCover
        active = surface.tilefrac .>= params.tthresh
        for (name, field, lo, hi) in (
                (:fveg, surface.fveg, Tf(0.02), Tf(0.99)),
                (:fves, surface.fves, Tf(0.02), Tf(0.99)),
                (:vfhp, surface.vfhp, Tf(0.02), Tf(0.99)),
                (:hcan, surface.hcan, Tf(1.0), Tf(100.0)),
                (:lai, surface.lai, Tf(0.05), Tf(10.0)),
            )
            bad = count(active .& ((field .< lo) .| (field .> hi)))
            bad == 0 || error("forest tile: $bad active cell(s) have $name outside [$lo, $hi]")
        end
    end

    return nothing
end

# Build the initial State from the finalized Surface.
function build_state(grid::Grid{Tf}, surface, params, substrate) where {Tf}
    GT = typeof(grid)
    st = State{GT, Matrix{Tf}, Matrix{Int}, Array{Tf, 3}}(; grid = grid)

    # Initialize soil states
    for k in 1:grid.Nsoil
        st.theta[k, :, :] .= params.fsat * surface.Vsat[:, :]
        st.Tsoil[k, :, :] .= params.Tprof
    end

    # Cap surface and soil temperatures on glacier ice
    if substrate isa IceSubstrate
        Tm = get_constants(Tf).Tm
        st.Tsrf .= min.(st.Tsrf, Tm)
        st.Tsoil .= min.(st.Tsoil, Tm)
    end

    return st
end

"""
    FSM(grid, landuse; arch = CPU(), params = Parameters, land_cover, substrate, ...schemes)
    FSM(Tf, landuse; ...)

Construct a ready-to-run FSM for one tile from a landuse domain. The tile is inferred from the
physics scheme types (`ForestCover`/`OpenCover`, `SoilSubstrate`/`IceSubstrate`). Each scheme
keyword accepts a scheme **type** (default-constructed at the grid precision) or a ready
**instance**. With a non-`CPU` `arch` the model is moved to that device.
"""
function FSM(
        grid::Grid{Tf}, landuse::Dict;
        arch::AbstractArchitecture = CPU(),
        params = Parameters{Tf}(),
        snow_albedo = PrognosticAlbedo{Tf}(grid),
        land_cover = OpenCover{Tf}(),
        substrate = SoilSubstrate{Tf}(),
        conductivity = DensityConductivity{Tf}(),
        fresh_snow_density = ElevationFreshSnowDensity{Tf}(),
        compaction = CrocusCompaction{Tf}(),
        hydrology = DensityBucketHydrology{Tf}(),
        layering = OriginalLayering{Tf}(),
        snow_fraction = PointSnowFraction{Tf}(),
    ) where {Tf}

    check_layer_thicknesses(grid)

    physics = (;
        snow_albedo = instantiate(snow_albedo, grid),
        land_cover = instantiate(land_cover, grid),
        substrate = instantiate(substrate, grid),
        conductivity = instantiate(conductivity, grid),
        fresh_snow_density = instantiate(fresh_snow_density, grid),
        compaction = instantiate(compaction, grid),
        hydrology = instantiate(hydrology, grid),
        layering = instantiate(layering, grid),
        snow_fraction = instantiate(snow_fraction, grid),
    )

    all(s -> s isa AbstractParameterization{Tf}, values(physics)) ||
        throw(ArgumentError("physics scheme precision does not match model Tf = $Tf"))
    for scheme in physics
        check_grid(scheme, grid.Nx, grid.Ny)
    end

    # Lock fresh snow density to a constant when the fixed scheme is selected
    if physics.fresh_snow_density isa FixedFreshSnowDensity
        params = reconstruct(params; rhof = params.rho0)
    end

    surface = build_surface(grid, landuse, physics.land_cover, physics.substrate, params)
    validate_surface(surface, params, physics.land_cover)
    state = build_state(grid, surface, params, physics.substrate)
    diag = Diagnostics{typeof(grid), Matrix{Tf}, Array{Tf, 3}}(; grid = grid)

    fsm = FSM(grid, params, surface, state, diag, physics)

    if !(arch isa CPU)
        fsm = on_architecture(arch, fsm)
    end
    return fsm
end

FSM(::Type{Tf}, landuse::Dict; kwargs...) where {Tf} = FSM(Grid(Tf, landuse), landuse; kwargs...)

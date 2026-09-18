# Test helper: build an FSM from an operational settings Dict (tile / physics / params) — the same
# Dict shape the FSM tests (and OSHDinternal) use.

using FlexibleSnowModelOSHD
import FlexibleSnowModelOSHD: Parameters, Surface, State, Diagnostics, check_layer_thicknesses

# Build a bare FSM (default/empty arrays, no landuse) for physics/transport/architecture unit tests
# that populate the state by hand.
function bare_fsm(
        grid::Grid{Tf};
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
    GT = typeof(grid)
    params = Parameters{Tf}()
    surface = Surface{GT, Matrix{Tf}, Matrix{Float64}, Matrix{Bool}}(; grid = grid)
    state = State{GT, Matrix{Tf}, Matrix{Int}, Array{Tf, 3}}(; grid = grid)
    diag = Diagnostics{GT, Matrix{Tf}, Array{Tf, 3}}(; grid = grid)
    physics = (;
        snow_albedo, land_cover, substrate, conductivity, fresh_snow_density,
        compaction, hydrology, layering, snow_fraction,
    )
    return FSM(grid, params, surface, state, diag, physics)
end

function build_fsm(arch::AbstractArchitecture, grid::Grid, landuse::Dict, settings::AbstractDict)
    tile = settings["tile"]
    tile in ("open", "forest", "glacier") ||
        error("tile requires open, forest or glacier (got tile = $tile)")

    physics = get(settings, "physics", Dict())
    land_cover = get(physics, "land_cover", tile == "forest" ? ForestCover : OpenCover)
    substrate = get(physics, "substrate", tile == "glacier" ? IceSubstrate : SoilSubstrate)
    schemes = (; (Symbol(k) => v for (k, v) in physics if !(k in ("land_cover", "substrate")))...)

    # Which cells to run: open runs everywhere; forest/glacier where the tile fraction >= tthresh
    # (0.1 is the test suite's threshold; reproduces the old tilefrac gate).
    Tf = eltype(grid)
    tthresh = Tf(0.1)
    active = tile == "open" ? trues(grid.Nx, grid.Ny) : (Tf.(landuse[tile]["data"]) .>= tthresh)

    fsm = FSM(grid, landuse; arch, active, land_cover, substrate, schemes...)

    # Per-cell Surface overrides (e.g. z0_snow). Scalar Parameters overrides are not used in the
    # FSM test suite, so they are intentionally unsupported here.
    for (k, v) in get(settings, "params", Dict())
        arr = getfield(fsm.surface, Symbol(k))
        v isa AbstractArray ? (arr .= eltype(arr).(v)) : fill!(arr, eltype(arr)(v))
    end
    return fsm
end

build_fsm(grid::Grid, landuse::Dict, settings::AbstractDict) = build_fsm(CPU(), grid, landuse, settings)

"""
$(TYPEDEF)

Static per-cell surface, terrain, canopy and soil properties, set at construction from the
landuse domain and read-only during the time loop.

# Fields

$(TYPEDFIELDS)
"""
@kwdef struct Surface{GT, MF, MF64, MB}
    grid::GT

    # Terrain
    "Grid elevation (m)"
    dem::MF = fill(NaN, grid.Nx, grid.Ny)
    "Grid cell size (m)"
    Ld::MF = fill(NaN, grid.Nx, grid.Ny)
    "Slope parameter (-)"
    slopemu::MF = fill(NaN, grid.Nx, grid.Ny)
    "Terrain correlation length (m)"
    xi::MF = fill(NaN, grid.Nx, grid.Ny)
    "Sky view fraction terrain (-)"
    fsky_terr::MF = fill(NaN, grid.Nx, grid.Ny)
    "Per-cell mask: run physics where true"
    active::MB = trues(grid.Nx, grid.Ny)

    # Snow and ground surface
    "Roughness length of snow (m)"
    z0_snow::MF = 0.002 * ones(grid.Nx, grid.Ny)
    "Snow-free roughness length (m)"
    z0sf::MF = 0.2 * ones(grid.Nx, grid.Ny)
    "Snow-free ground albedo (-)"
    alb0::MF = 0.2 * ones(grid.Nx, grid.Ny)

    # Forest
    "Vegetation area index (-)"
    VAI::MF = fill(NaN, grid.Nx, grid.Ny)
    "Leaf area index (-)"
    lai::MF = fill(NaN, grid.Nx, grid.Ny)
    "Canopy cover fraction (-)"
    fveg::MF = fill(NaN, grid.Nx, grid.Ny)
    "Stand-scale canopy cover fraction (-)"
    fves::MF = fill(NaN, grid.Nx, grid.Ny)
    "Sky view fraction (-)"
    fsky::MF = fill(NaN, grid.Nx, grid.Ny)
    "Hemispherical sky-view fraction incl. canopy (-)"
    vfhp::MF = fill(NaN, grid.Nx, grid.Ny)
    "Canopy height (m)"
    hcan::MF = fill(NaN, grid.Nx, grid.Ny)
    "Canopy heat capacity (J/K/m^2)"
    canh::MF = fill(NaN, grid.Nx, grid.Ny)
    "Canopy snow capacity (kg/m^2)"
    scap::MF = fill(NaN, grid.Nx, grid.Ny)
    "Canopy transmissivity (-)"
    trcn::MF = ones(grid.Nx, grid.Ny)

    # Precipitation
    "Precipitation multiplier reverting open-area correction (-)"
    pmultf::MF = fill(NaN, grid.Nx, grid.Ny)
    "Precipitation multiplier (-); Float64 legacy"
    prec_multi::MF64 = fill(NaN, grid.Nx, grid.Ny)

    # Soil
    "Soil clay fraction (-)"
    fcly::MF = 0.3 * ones(grid.Nx, grid.Ny)
    "Soil sand fraction (-)"
    fsnd::MF = 0.6 * ones(grid.Nx, grid.Ny)
    "Clapp-Hornberger exponent (-)"
    b::MF = zeros(grid.Nx, grid.Ny)
    "Volumetric heat capacity of dry soil (J/K/m^3)"
    hcap_soil::MF = zeros(grid.Nx, grid.Ny)
    "Thermal conductivity of dry soil (W/m/K)"
    hcon_soil::MF = zeros(grid.Nx, grid.Ny)
    "Saturated soil water pressure (m)"
    sathh::MF = zeros(grid.Nx, grid.Ny)
    "Volumetric soil moisture at saturation (-)"
    Vsat::MF = zeros(grid.Nx, grid.Ny)
    "Volumetric soil moisture at critical point (-)"
    Vcrit::MF = zeros(grid.Nx, grid.Ny)
end

@adapt_structure Surface

function Base.show(io::IO, x::Surface)
    grid = x.grid
    print(io, nameof(typeof(x)), '\n')
    print(io, "├── Precision: ", eltype(grid), '\n')
    print(io, "├── Grid: ", grid.Nx, " × ", grid.Ny, '\n')
    print(io, "└── Fields: ", length(propertynames(x)) - 1, '\n')
    return nothing
end

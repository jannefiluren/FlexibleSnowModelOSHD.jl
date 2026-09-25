"""
$(TYPEDEF)

Prognostic model state carried across time steps: the snow, canopy and soil variables
updated in place by `step!`.

# Fields

$(TYPEDFIELDS)
"""
@kwdef struct State{GT, MF, MI, AF}
    grid::GT
    "Snow albedo (-)"
    albs::MF = 0.85 * ones(grid.Nx, grid.Ny)
    "Snow layer thicknesses (m)"
    Ds::AF = zeros(grid.Nsmax, grid.Nx, grid.Ny)
    "Number of snow layers"
    Nsnow::MI = zeros(Int, grid.Nx, grid.Ny)
    "Canopy air space humidity (kg/kg)"
    Qcan::MF = zeros(grid.Nx, grid.Ny)
    "Ice content of snow layers (kg/m^2)"
    Sice::AF = zeros(grid.Nsmax, grid.Nx, grid.Ny)
    "Liquid content of snow layers (kg/m^2)"
    Sliq::AF = zeros(grid.Nsmax, grid.Nx, grid.Ny)
    "Snow mass on vegetation (kg/m^2)"
    Sveg::MF = zeros(grid.Nx, grid.Ny)
    "Canopy air space temperature (K)"
    Tcan::MF = 285 * ones(grid.Nx, grid.Ny)
    "Volumetric moisture content of soil layers (-)"
    theta::AF = zeros(grid.Nsoil, grid.Nx, grid.Ny)
    "Snow layer temperatures (K)"
    Tsnow::AF = 273.15 * ones(grid.Nsmax, grid.Nx, grid.Ny)
    "Soil layer temperatures (K)"
    Tsoil::AF = 285 * ones(grid.Nsoil, grid.Nx, grid.Ny)
    "Surface skin temperature (K)"
    Tsrf::MF = 285 * ones(grid.Nx, grid.Ny)
    "Snow cover fraction (-)"
    fsnow::MF = zeros(grid.Nx, grid.Ny)
    "Vegetation temperature (K)"
    Tveg::MF = 285 * ones(grid.Nx, grid.Ny)
    "Min snow depth at time of swemin (m)"
    snowdepthmin::MF = zeros(grid.Nx, grid.Ny)
    "Max snow depth at time of swemax (m)"
    snowdepthmax::MF = zeros(grid.Nx, grid.Ny)
    "Snow depth over last 14 days (m)"
    snowdepthhist::AF = zeros(14, grid.Nx, grid.Ny)
    "Minimum SWE during the season (kg/m^2)"
    swemin::MF = zeros(grid.Nx, grid.Ny)
    "Maximum SWE during the season (kg/m^2)"
    swemax::MF = zeros(grid.Nx, grid.Ny)
    "SWE over last 14 days (kg/m^2)"
    swehist::AF = zeros(14, grid.Nx, grid.Ny)
    "Historical past wetting of a layer (-)"
    histowet::AF = zeros(grid.Nsmax, grid.Nx, grid.Ny)
end

@adapt_structure State

function Base.show(io::IO, x::State)
    grid = x.grid
    print(io, nameof(typeof(x)), '\n')
    print(io, "├── Precision: ", eltype(grid), '\n')
    print(io, "├── Grid: ", grid.Nx, " × ", grid.Ny, '\n')
    print(io, "└── Fields: ", length(propertynames(x)) - 1, '\n')
    return nothing
end

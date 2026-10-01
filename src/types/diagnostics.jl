"""
$(TYPEDEF)

Per-cell diagnostic fields recomputed within a time step: intermediate fluxes and
properties passed between the physics processes. Grouped by the process that fills them.

# Fields

$(TYPEDFIELDS)
"""
@kwdef struct Diagnostics{GT, MF, AF}
    grid::GT

    # drive
    "Saturation vapour pressure (Pa)"
    es::MF = zeros(grid.Nx, grid.Ny)
    "Specific humidity (kg/kg)"
    Qa::MF = zeros(grid.Nx, grid.Ny)
    "Wind speed with lower bound applied (m/s)"
    Uaeff::MF = zeros(grid.Nx, grid.Ny)
    "Snowfall reaching the surface (kg/m^2/s)"
    Sfeff::MF = zeros(grid.Nx, grid.Ny)

    # radiation
    "Net shortwave absorbed by vegetation (W/m^2)"
    SWveg::MF = zeros(grid.Nx, grid.Ny)
    "Net shortwave absorbed by the surface (W/m^2)"
    SWsrf::MF = zeros(grid.Nx, grid.Ny)
    "Subcanopy incoming shortwave (W/m^2)"
    SWsci::MF = zeros(grid.Nx, grid.Ny)
    "Incoming longwave used in the energy balance (W/m^2)"
    LWeff::MF = zeros(grid.Nx, grid.Ny)

    # thermal
    "Thermal conductivity of snow (W/m/K)"
    ksnow::AF = zeros(grid.Nsmax, grid.Nx, grid.Ny)
    "Areal heat capacity of soil (J/K/m^2)"
    csoil::AF = zeros(grid.Nsoil, grid.Nx, grid.Ny)
    "Thermal conductivity of soil (W/m/K)"
    ksoil::AF = zeros(grid.Nsoil, grid.Nx, grid.Ny)
    "Surface moisture conductance (m/s)"
    gs1::MF = zeros(grid.Nx, grid.Ny)
    "Surface layer thickness (m)"
    Ds1::MF = zeros(grid.Nx, grid.Ny)
    "Surface layer temperature (K)"
    Ts1::MF = zeros(grid.Nx, grid.Ny)
    "Surface thermal conductivity (W/m/K)"
    ks1::MF = zeros(grid.Nx, grid.Ny)
    "Vegetation temperature at start of timestep (K)"
    Tveg0::MF = zeros(grid.Nx, grid.Ny)

    # surface_exchange_coefficients
    "Eddy diffusivity for heat to the atmosphere (m/s)"
    KH::MF = zeros(grid.Nx, grid.Ny)
    "Eddy diffusivity from the canopy air space (m/s)"
    KHa::MF = zeros(grid.Nx, grid.Ny)
    "Eddy diffusivity for heat from the ground (m/s)"
    KHg::MF = zeros(grid.Nx, grid.Ny)
    "Eddy diffusivity for heat from vegetation (m/s)"
    KHv::MF = zeros(grid.Nx, grid.Ny)
    "Eddy diffusivity for water from the ground (m/s)"
    KWg::MF = zeros(grid.Nx, grid.Ny)
    "Eddy diffusivity for water from vegetation (m/s)"
    KWv::MF = zeros(grid.Nx, grid.Ny)
    "Wind speed in canopy layer (m/s)"
    Usc::MF = zeros(grid.Nx, grid.Ny)

    # surface_energy_balance
    "Moisture flux from the surface (kg/m^2/s)"
    Esrf::MF = zeros(grid.Nx, grid.Ny)
    "Moisture flux from vegetation (kg/m^2/s)"
    Eveg::MF = zeros(grid.Nx, grid.Ny)
    "Heat flux into the surface (W/m^2)"
    G::MF = zeros(grid.Nx, grid.Ny)
    "Sensible heat flux to the atmosphere (W/m^2)"
    H::MF = zeros(grid.Nx, grid.Ny)
    "Sensible heat flux from the surface (W/m^2)"
    Hsrf::MF = zeros(grid.Nx, grid.Ny)
    "Latent heat flux to the atmosphere (W/m^2)"
    LE::MF = zeros(grid.Nx, grid.Ny)
    "Latent heat flux from the surface (W/m^2)"
    LEsrf::MF = zeros(grid.Nx, grid.Ny)
    "Subcanopy incoming longwave (W/m^2)"
    LWsci::MF = zeros(grid.Nx, grid.Ny)
    "Net longwave absorbed by vegetation (W/m^2)"
    LWveg::MF = zeros(grid.Nx, grid.Ny)
    "Surface melt rate (kg/m^2/s)"
    Melt::MF = zeros(grid.Nx, grid.Ny)
    "Net radiation (W/m^2)"
    Rnet::MF = zeros(grid.Nx, grid.Ny)
    "Net radiation at surface (W/m^2)"
    Rsrf::MF = zeros(grid.Nx, grid.Ny)

    # canopy
    "Canopy interception (kg/m^2)"
    intcpt::MF = zeros(grid.Nx, grid.Ny)
    "Sublimation from vegetation (kg/m^2)"
    Sbveg::MF = zeros(grid.Nx, grid.Ny)
    "Snow mass unloaded from canopy (kg/m^2)"
    unload::MF = zeros(grid.Nx, grid.Ny)

    # snow
    "Heat flux into soil (W/m^2)"
    Gsoil::MF = zeros(grid.Nx, grid.Ny)
    "Total runoff (kg/m^2)"
    Roff::MF = zeros(grid.Nx, grid.Ny)
    "Runoff from snowmelt at base of snow (kg/m^2)"
    meltflux_out::MF = zeros(grid.Nx, grid.Ny)
    "Sublimation from the snow surface (kg/m^2)"
    Sbsrf::MF = zeros(grid.Nx, grid.Ny)
    "Bare soil runoff (kg/m^2)"
    Roff_bare::MF = zeros(grid.Nx, grid.Ny)
    "Runoff at base of snow (kg/m^2)"
    Roff_snow::MF = zeros(grid.Nx, grid.Ny)
    "Snow depth at start of timestep (m)"
    snowdepth0::MF = zeros(grid.Nx, grid.Ny)
    "Ice content at start of timestep (kg/m^2)"
    Sice0::MF = zeros(grid.Nx, grid.Ny)

    # snow_layering
    "Snow layer thickness at start of timestep (m)"
    Ds0::MF = zeros(grid.Nx, grid.Ny)
end

@adapt_structure Diagnostics

function Base.show(io::IO, x::Diagnostics)
    grid = x.grid
    print(io, nameof(typeof(x)), '\n')
    print(io, "├── Precision: ", eltype(grid), '\n')
    print(io, "├── Grid: ", grid.Nx, " × ", grid.Ny, '\n')
    print(io, "└── Fields: ", length(propertynames(x)) - 1, '\n')
    return nothing
end

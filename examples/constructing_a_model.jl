# # Constructing a model
#
# This example builds an `FSM` model from scratch on a small grid, runs it with a single
# time step of meteorological forcing, and reads out the snow depth — first over open
# terrain, then with a forest canopy.

using FlexibleSnowModelOSHD
using Dates

# ## Meteorological forcing
#
# Create the forcing for one time step on a `2 × 3` grid at `Float32` precision.

Nx, Ny = 2, 3
precision = Float32

met = MET{precision}(
    Sdir = fill(200, Nx, Ny),         # Direct shortwave radiation per inclined surface area (W/m^2)
    Sdird = fill(200, Nx, Ny),        # Direct shortwave radiation per horizontal surface area (W/m^2)
    Sdif = fill(100, Nx, Ny),         # Diffuse shortwave radiation (W/m^2)
    LW = fill(300, Nx, Ny),           # Incoming longwave radiation (W/m^2)
    Sf = fill(0.01, Nx, Ny),          # Snowfall rate (kg/m^2/s)
    Rf = fill(0, Nx, Ny),             # Rainfall rate (kg/m^2/s)
    Sf24h = fill(10, Nx, Ny),         # Total snowfall over 24h (kg/m^2)
    Ta = fill(270, Nx, Ny),           # Air temperature (K)
    RH = fill(80, Nx, Ny),            # Relative humidity (%)
    Ua = fill(3, Nx, Ny),             # Wind speed (m/s)
    Ps = fill(100000, Nx, Ny),        # Surface air pressure (Pa)
    Tv = fill(0.5, Nx, Ny),           # Time-varying transmissivity for direct shortwave radiation (-)
)

# ## Landuse
#
# The static per-cell terrain, canopy and precipitation properties the model needs.

landuse = Dict(
    "skyvf" => Dict("data" => fill(0.95, Nx, Ny)),
    "elevation" => Dict("data" => fill(2540.0, Nx, Ny)),
    "prec_multi" => Dict("data" => ones(Float64, Nx, Ny)),
    "slopemu" => Dict("data" => ones(Float64, Nx, Ny)),
    "xi" => Dict("data" => fill(150.0, Nx, Ny)),
    "Ld" => Dict("data" => fill(250.0, Nx, Ny)),
    "forest" => Dict("data" => fill(0.6, Nx, Ny)),
    "glacier" => Dict("data" => fill(0.5, Nx, Ny)),
    "fveg" => Dict("data" => fill(0.6, Nx, Ny)),
    "fves" => Dict("data" => fill(0.6, Nx, Ny)),
    "hcan" => Dict("data" => fill(20.0, Nx, Ny)),
    "lai" => Dict("data" => fill(2.5, Nx, Ny)),
    "vfhp" => Dict("data" => fill(0.5, Nx, Ny)),
)

# ## Open terrain
#
# Build the model with default settings, run one time step, and compute total snow depth.

fsm = FSM(precision, landuse)
step!(fsm, met, DateTime(2023, 12, 1, 12))
snowdepth = dropdims(sum(fsm.state.Ds, dims = 1), dims = 1)

# ## Forest cover
#
# The same domain with a `ForestCover` land cover, so the canopy modifies the surface forcing.

land_cover = ForestCover{precision}()
fsm = FSM(precision, landuse, land_cover = land_cover)
step!(fsm, met, DateTime(2023, 12, 1, 12))
snowdepth = dropdims(sum(fsm.state.Ds, dims = 1), dims = 1)

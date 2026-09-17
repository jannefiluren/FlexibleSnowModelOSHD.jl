# FlexibleSnowModelOSHD.jl

[ci-img]: https://github.com/jannefiluren/FlexibleSnowModelOSHD.jl/actions/workflows/CI.yml/badge.svg?branch=main
[ci-url]: https://github.com/jannefiluren/FlexibleSnowModelOSHD.jl/actions/workflows/CI.yml?query=branch%3Amain

[codecov-img]: https://codecov.io/gh/jannefiluren/FlexibleSnowModelOSHD.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/jannefiluren/FlexibleSnowModelOSHD.jl

[runic-img]: https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-black
[runic-url]: https://github.com/fredrikekre/Runic.jl

[zenodo-img]: https://img.shields.io/badge/DOI-10.5281/zenodo.21322824-blue.svg
[zenodo-url]: https://doi.org/10.5281/zenodo.21322824

 [![][ci-img]][ci-url] [![][codecov-img]][codecov-url] [![][runic-img]][runic-url] [![DOI][zenodo-img]][zenodo-url]

A Julia implementation of the **Flexible Snow Model (FSM)** for the **Operational Snow Hydrological Service (OSHD)** at SLF. This package provides a comprehensive snow physics model for simulating snow accumulation and melt processes in complex terrain.

## Overview

FlexibleSnowModelOSHD is a multi-layer snow model that simulates:

- **Snow accumulation and ablation** with detailed physics-based processes
- **Multi-layer snow structure** with dynamic layer evolution 
- **Energy balance** including shortwave/longwave radiation, turbulent fluxes, and ground heat transfer
- **Snow hydraulics** with configurable drainage schemes and liquid water retention
- **Forest canopy interactions** including snow interception, unloading, and subcanopy processes
- **Fractional snow cover** using multiple parameterization approaches

The model is designed for operational snow forecasting applications and supports both point-scale and distributed (gridded) simulations across various surface types including open areas, forests, and glaciers.

## Installation

### Prerequisites
- Julia 1.10 or higher
- Required packages are specified in `Project.toml`

### Installation Steps

1. **Clone the repository:**
   ```bash
   git clone https://github.com/jannefiluren/FlexibleSnowModelOSHD.jl
   cd FlexibleSnowModelOSHD
   ```

2. **Activate the package environment:**
   ```julia
   using Pkg
   Pkg.activate(".")
   Pkg.instantiate()
   ```

3. **Load the package:**
   ```julia
   using FlexibleSnowModelOSHD
   ```

## Examples

A simulation representing an open site can be run from the terminal by:

```julia
include("script/run_open_station_example.jl")
```

while a corresponding simulation for a forested site can be run by:

```julia
include("script/run_forest_station_example.jl")
```

## Constructing a model

Load packages:

```julia
using FlexibleSnowModelOSHD
using Dates
```

Create meteorological forcing data for one time step:

```julia
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
```

Create landuse data needed for the model:

```julia
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
```

Setup the model using default settings, run the model one time step and compute total snow depth:

```julia
fsm = FSM(precision, landuse)
step!(fsm, met, DateTime(2023, 12, 1, 12))
snowdepth = dropdims(sum(fsm.state.Ds, dims=1), dims=1)
```

Setup the model for domain with forest cover:

```julia
land_cover = ForestCover{precision}()
fsm = FSM(precision, landuse, land_cover = land_cover)
step!(fsm, met, DateTime(2023, 12, 1, 12))
snowdepth = dropdims(sum(fsm.state.Ds, dims=1), dims=1)
```


## Package Structure

```
FlexibleSnowModelOSHD.jl/
├── src/                              # Source code
│   ├── FlexibleSnowModelOSHD.jl      # Main module
│   ├── parameters.jl                 # Physical constants
│   ├── types.jl                      # Model data structures
│   ├── architectures.jl              # CPU/GPU architecture abstraction
│   ├── scheme_construction.jl        # Scheme construction helpers
│   ├── construct.jl                  # Model construction from a landuse domain
│   ├── step.jl                       # Main physics time step
│   ├── numerics/                     # Numerical utilities
│   ├── processes/                    # Kernel functions called from step
│   ├── parameterizations/            # Swappable physics schemes
│   │   ├── albedo.jl                 # Snow albedo schemes
│   │   ├── compaction.jl             # Snow compaction schemes
│   │   ├── conductivity.jl           # Snow thermal conductivity schemes
│   │   ├── fresh_snow_density.jl     # Fresh snow density schemes
│   │   ├── hydrology.jl              # Snow liquid-water (hydrology) schemes
│   │   ├── land_cover.jl             # Open/forest land cover and canopy schemes
│   │   ├── layering.jl               # Snow layering schemes
│   │   ├── snow_cover_fraction.jl    # Snow cover fraction schemes
│   │   ├── stability.jl              # Atmospheric stability correction schemes
│   │   └── substrate.jl              # Soil / glacier-ice substrate schemes
│   └── transport/                    # Snow transport routines
├── script/                           # Simulation scripts
└── test/                             # Unit tests and regression tests
```

# FlexibleSnowModelOSHD.jl

[ci-img]: https://github.com/jannefiluren/FlexibleSnowModelOSHD.jl/actions/workflows/CI.yml/badge.svg?branch=main
[ci-url]: https://github.com/jannefiluren/FlexibleSnowModelOSHD.jl/actions/workflows/CI.yml?query=branch%3Amain

[codecov-img]: https://codecov.io/gh/jannefiluren/FlexibleSnowModelOSHD.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/jannefiluren/FlexibleSnowModelOSHD.jl

[runic-img]: https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-black
[runic-url]: https://github.com/fredrikekre/Runic.jl

[zenodo-img]: https://img.shields.io/badge/DOI-10.5281/zenodo.21322824-blue.svg
[zenodo-url]: https://doi.org/10.5281/zenodo.21322824

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://jannefiluren.github.io/FlexibleSnowModelOSHD.jl/dev/

 [![][docs-dev-img]][docs-dev-url] [![][ci-img]][ci-url] [![][codecov-img]][codecov-url] [![][runic-img]][runic-url] [![DOI][zenodo-img]][zenodo-url]

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

## Package Structure

```
FlexibleSnowModelOSHD.jl/
├── src/                              # Source code
│   ├── FlexibleSnowModelOSHD.jl      # Main module
│   ├── parameters.jl                 # Physical constants
│   ├── types/                        # Core model data structures
│   │   ├── grid.jl                   # Grid definition
│   │   ├── parameters.jl             # Scalar model parameters
│   │   ├── surface.jl                # Static surface / terrain / canopy / soil properties
│   │   ├── state.jl                  # Prognostic model state
│   │   ├── diagnostics.jl            # Model diagnostic fields
│   │   ├── met.jl                    # Meteorological forcing
│   │   └── fsm.jl                    # FSM container
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
├── examples/                         # Runnable example scripts (rendered into the docs)
└── test/                             # Unit tests and regression tests
```

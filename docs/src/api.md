# API reference

Documentation grouped by topic. Each section displays docstrings from the source files, in
an order that follows the model structure.

## Core types

```@autodocs
Modules = [FlexibleSnowModelOSHD]
Pages = [
    "types/fsm.jl",
    "types/grid.jl",
    "types/parameters.jl",
    "types/surface.jl",
    "types/state.jl",
    "types/diagnostics.jl",
    "types/met.jl",
]
```

## Model construction and stepping

```@autodocs
Modules = [FlexibleSnowModelOSHD]
Pages = ["construct.jl", "step.jl", "scheme_construction.jl"]
```

## Physics schemes

```@autodocs
Modules = [FlexibleSnowModelOSHD]
Pages = [
    "parameterizations/albedo.jl",
    "parameterizations/conductivity.jl",
    "parameterizations/compaction.jl",
    "parameterizations/fresh_snow_density.jl",
    "parameterizations/hydrology.jl",
    "parameterizations/layering.jl",
    "parameterizations/snow_cover_fraction.jl",
    "parameterizations/stability.jl",
    "parameterizations/land_cover.jl",
    "parameterizations/substrate.jl",
]
```

## Physical processes

```@autodocs
Modules = [FlexibleSnowModelOSHD]
Pages = [
    "processes/drive.jl",
    "processes/radiation.jl",
    "processes/surface_energy_balance.jl",
    "processes/surface_exchange_coefficients.jl",
    "processes/thermal.jl",
    "processes/snow.jl",
    "processes/soil.jl",
    "processes/canopy_mass_balance.jl",
]
```

## Numerics

```@autodocs
Modules = [FlexibleSnowModelOSHD]
Pages = ["numerics/qsat.jl", "numerics/tridiag.jl", "numerics/ludcmp.jl", "numerics/reductions.jl"]
```

## Architectures (CPU / GPU)

```@autodocs
Modules = [FlexibleSnowModelOSHD]
Pages = ["architectures.jl"]
```

## Snow transport

```@autodocs
Modules = [FlexibleSnowModelOSHD]
Pages = [
    "transport/transport.jl",
    "transport/transport_setup.jl",
    "transport/snowtran3d.jl",
    "transport/snowtran3d_julia.jl",
    "transport/snowslide.jl",
    "transport/snowslide_julia.jl",
    "transport/relayer.jl",
]
```

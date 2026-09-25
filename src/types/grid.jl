"""
$(TYPEDEF)

Grid definition with a given floating point precision, size and snow/soil layer
thicknesses. Construct with `Grid(Tf; kwargs...)`:

```jldoctest
using FlexibleSnowModelOSHD

grid = Grid(Float32, Nx = 10, Ny = 5)

# output
Grid
├── Precision: Float32
├── Nx: 10
├── Ny: 5
├── Dzsnow: [0.1, 0.2, 0.4]
└── Dzsoil: [0.1, 0.2, 0.4, 0.8]
```

# Fields

$(TYPEDFIELDS)
"""
@kwdef struct Grid{Tf, VF <: AbstractVector{Tf}}
    "Maximum snow layer thicknesses (m)"
    Dzsnow::VF = Tf[0.1, 0.2, 0.4]
    "Maximum soil layer thicknesses (m)"
    Dzsoil::VF = Tf[0.1, 0.2, 0.4, 0.8]
    "Number of snow layers"
    Nsmax::Int = length(Dzsnow)
    "Number of soil layers"
    Nsoil::Int = length(Dzsoil)
    "First array dimension (rows)"
    Nx::Int = 1
    "Second array dimension (columns)"
    Ny::Int = 1
end

function Grid(::Type{Tf}; kwargs...) where {Tf}
    grid = Grid{Tf, Vector{Tf}}(; kwargs...)
    check_layer_thicknesses(grid)
    return grid
end

Base.eltype(::Grid{Tf}) where {Tf} = Tf
Base.eltype(::Type{<:Grid{Tf}}) where {Tf} = Tf

"""
$(TYPEDSIGNATURES)

Verify that the first snow layer can grow at least as thick as the top soil layer.
The surface layer in `thermal!` blends snow and soil over a depth of `Dzsoil[1]`, so a
thinner first snow layer leaves `Ts1` and `ks1` contaminated by soil however deep the
snowpack gets.
"""
function check_layer_thicknesses(grid::Grid)
    grid.Dzsnow[1] >= grid.Dzsoil[1] || throw(
        ArgumentError(
            "Dzsnow[1] = $(grid.Dzsnow[1]) must be at least Dzsoil[1] = $(grid.Dzsoil[1])"
        )
    )
    return nothing
end

@adapt_structure Grid

function Base.show(io::IO, grid::Grid{Tf}) where {Tf}
    print(io, "Grid", '\n')
    print(io, "├── Precision: ", Tf, '\n')
    print(io, "├── Nx: ", grid.Nx, '\n')
    print(io, "├── Ny: ", grid.Ny, '\n')
    print(io, "├── Dzsnow: ", "[" * join(grid.Dzsnow, ", ") * "]", '\n')
    print(io, "└── Dzsoil: ", "[" * join(grid.Dzsoil, ", ") * "]", '\n')
    return nothing
end

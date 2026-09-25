# Helpers for materializing physics parameterizations against a grid, plus the immutable-struct
# update used to rebuild `Parameters`. Used by the FSM constructor (construct.jl) and by callers
# that translate an operational config into constructor arguments.

"""
$(TYPEDSIGNATURES)

Materialize a parameterization parameter as an `Nx` by `Ny` array of element type `Tf`.
A scalar is broadcast over the whole grid; an array is converted element-wise.
"""
grid_array(::Type{Tf}, x::Number, Nx, Ny) where {Tf} = fill(Tf(x), Nx, Ny)
grid_array(::Type{Tf}, x::AbstractArray, Nx, Ny) where {Tf} = convert(Array{Tf, 2}, x)

"""
$(TYPEDSIGNATURES)

Assert that any grid-shaped parameter held by `scheme` matches the `Nx` by `Ny` model grid.
The fallback accepts anything, so a parameterization built only from scalars needs no method.

Without this a mismatch is not caught at setup; it surfaces later as a `BoundsError` from
inside a kernel, which says nothing about the actual cause.
"""
check_grid(scheme, Nx, Ny) = nothing

function check_grid(scheme::AbstractParameterization, Nx, Ny)
    for name in fieldnames(typeof(scheme))
        value = getfield(scheme, name)
        value isa AbstractArray || continue
        size(value) == (Nx, Ny) || throw(
            DimensionMismatch(
                "$(nameof(typeof(scheme))) field `$name` is $(size(value)) but the model grid is ($Nx, $Ny)"
            )
        )
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Return a physics parameterization ready for the model: a scheme **type** is default-constructed at
the grid's precision (`Scheme{eltype(grid)}(grid)`), while a ready-made **instance** is returned
unchanged. Non-default schemes are constructed by the caller, e.g.
`PrognosticAlbedo{Float32}(grid; adm = 200, adc = my_array)`.
"""
instantiate(scheme::Type, grid) = scheme{eltype(grid)}(grid)
instantiate(scheme, grid) = scheme

"""
$(TYPEDSIGNATURES)

Copy the immutable struct `x` with the named fields replaced. `Parameters` is rebuilt
rather than mutated so that it stays isbits and can cross into a kernel by value.
"""
function reconstruct(x::T; kwargs...) where {T}
    names = fieldnames(T)
    fields = NamedTuple{names}(map(f -> getfield(x, f), names))
    return T(; merge(fields, NamedTuple(kwargs))...)
end

function Base.show(io::IO, p::AbstractParameterization; indent = "")
    print(io, nameof(typeof(p)), '\n')
    fields = filter(f -> !(getfield(p, f) isa Grid), propertynames(p))
    for (i, field) in enumerate(fields)
        last = i == length(fields)
        start = indent * (last ? "└── " : "├── ")
        values = getfield(p, field)
        if values isa AbstractParameterization
            print(io, start, field, ": ")
            show(io, values; indent = indent * (last ? "    " : "│   "))
        elseif values isa AbstractArray && ndims(values) > 0
            print(io, start, field, ": ", extrema(values), '\n')
        else
            print(io, start, field, ": ", values, '\n')
        end
    end
    return nothing
end

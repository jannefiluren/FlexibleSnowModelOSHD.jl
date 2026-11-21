# One-time setup for the snow-transport operators: build a SnowTransport workspace from the
# model grid + a landuse dictionary, writing into the standalone workspace rather than onto FSM.
# Both the Julia and Fortran transport paths consume this same workspace, so any array computed
# here is shared and cannot make the two paths disagree.

"""
$(TYPEDSIGNATURES)

Allocate and initialise a [`SnowTransport`](@ref) workspace for `fsm`:

- loads `slope` (deg) and the `forest`/`glacier` fractions from the `landuse` dictionary,
- derives the snow holding depth `Shd` from slope,
- sorts the DEM indices (highest → lowest) for the SnowSlide processing order.

The vegetation snow-holding depth (0.1 m) and the tuning constants keep their
[`SnowTransport`](@ref) defaults, which match `deps/MODULES.F90`.

`wind`/`slide` select which processes [`transport!`](@ref) runs; `use_fortran` chooses the
Fortran ccall path over the Julia port. Transport runs in the open tile only, so the
workspace is built with the open-tile fractions: `tilefrac = 1 − forest − glacier`, with
glaciers counted as open (`openfrac = tilefrac + glacierfrac`). Transport is CPU-only, so
`fsm` must hold host arrays.
"""
function setup_transport(
        fsm::FSM{Tf}, landuse::Dict;
        wind::Bool = false, slide::Bool = false, use_fortran::Bool = false,
    ) where {Tf}

    Nx = Int(fsm.grid.Nx)
    Ny = Int(fsm.grid.Ny)

    w = SnowTransport{Tf}(Nx = Nx, Ny = Ny)
    w.wind = wind
    w.slide = slide
    w.use_fortran = use_fortran

    # Open fractions; transport weights mass by openfrac for the tile recombination
    forestfrac = haskey(landuse, "forest") ? Tf.(landuse["forest"]["data"]) : zeros(Tf, Nx, Ny)
    if haskey(landuse, "glacier")
        w.glacierfrac .= Tf.(landuse["glacier"]["data"])
    end
    w.tilefrac .= Tf(1) .- forestfrac .- w.glacierfrac
    w.openfrac .= w.tilefrac .+ w.glacierfrac

    # Slope (deg) — required
    haskey(landuse, "slope") ||
        throw(ArgumentError("setup_transport requires 'slope' data in the landuse dictionary"))
    w.slope .= Tf.(landuse["slope"]["data"])

    # Snow holding depth from slope: normal-to-slope value projected to the vertical.
    # Slope floored at 10 deg to avoid the divergence of the power law near 0.
    slope_thres = max.(w.slope, Tf(10))
    shd_norm = Tf(3178.4) .* slope_thres .^ Tf(-1.998)
    w.Shd .= shd_norm .* max.(cosd.(slope_thres), Tf(0.001))

    # Processing order for SnowSlide: cells from highest to lowest elevation
    sort_dem_indices!(w.index_sorted_dem, fsm.surface.dem)

    return w
end

"""
$(TYPEDSIGNATURES)

Sort DEM indices from highest to lowest elevation for the SnowSlide processing order. Fills
`index_sorted_dem`, an `(Nx*Ny, 2)` array whose rows are `(i, j)` pairs ordered by decreasing
`dem`.
"""
function sort_dem_indices!(index_sorted_dem::Matrix{Int}, dem::AbstractMatrix{Tf}) where {Tf <: Real}
    ind = sortperm(vec(dem), rev = true)

    rows, cols = size(dem)
    rows2D = repeat((1:rows), outer = (1, cols))
    cols2D = repeat((1:cols)', outer = (rows, 1))

    index_sorted_dem[:, 1] = Int.(rows2D[ind])
    index_sorted_dem[:, 2] = Int.(cols2D[ind])

    return nothing
end

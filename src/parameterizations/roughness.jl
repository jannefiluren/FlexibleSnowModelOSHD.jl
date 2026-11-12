"""
$(TYPEDEF)

Snow roughness length scheme using the per-cell `z0_snow` of the surface unchanged.

```jldoctest
using FlexibleSnowModelOSHD

FixedSnowRoughness{Float32}()

# output
FixedSnowRoughness
```
"""
struct FixedSnowRoughness{Tf} <: AbstractSnowRoughness{Tf} end

"""
$(TYPEDEF)

Snow roughness length scheme scaling the per-cell `z0_snow` with the liquid water content of the
top snow layer: from `z0_snow` for dry snow towards `z0_wet_ratio * z0_snow`, reached in full at the
volumetric liquid water content `theta_z0_wet`.

```jldoctest
using FlexibleSnowModelOSHD

WetSnowRoughness{Float32}(; z0_wet_ratio = 1.5)

# output
WetSnowRoughness
├── z0_wet_ratio: 1.5
└── theta_z0_wet: 0.01
```

# Fields

$(TYPEDFIELDS)
"""
@kwdef struct WetSnowRoughness{Tf} <: AbstractSnowRoughness{Tf}
    "Ratio of wet- to dry-snow roughness length; 1 disables the scaling (-)"
    z0_wet_ratio::Tf = 1
    "Volumetric liquid water content of the top snow layer at which `z0_wet_ratio` applies in full (-)"
    theta_z0_wet::Tf = 0.01
end

FixedSnowRoughness{Tf}(grid::Grid; kwargs...) where {Tf} = FixedSnowRoughness{Tf}()
WetSnowRoughness{Tf}(grid::Grid; kwargs...) where {Tf} = WetSnowRoughness{Tf}(; kwargs...)

"""
    snow_roughness(scheme, i, j, state, surface)

Roughness length of snow at cell `(i, j)`. Implemented for every `AbstractSnowRoughness`.
"""
function snow_roughness end

@inline snow_roughness(::FixedSnowRoughness, i, j, state, surface) = surface.z0_snow[i, j]

@inline function snow_roughness(scheme::WetSnowRoughness{Tf}, i, j, state, surface) where {Tf}
    @unpack_constants(Tf)
    (; Ds, Sliq, fsnow) = state
    (; z0_snow) = surface
    (; z0_wet_ratio, theta_z0_wet) = scheme

    if Ds[1, i, j] * fsnow[i, j] <= eps(Tf)
        return z0_snow[i, j]
    end
    # Ds is the depth over the snow-covered fraction, Sliq the cell mean
    theta = Sliq[1, i, j] / (rho_wat * Ds[1, i, j] * fsnow[i, j])
    return z0_snow[i, j] * (one(Tf) + (z0_wet_ratio - one(Tf)) * min(theta / theta_z0_wet, one(Tf)))
end

"""
    ground_roughness(snow_fraction, snow_roughness, i, j, state, surface)

Roughness length of the ground at cell `(i, j)`: the [`snow_roughness`](@ref) where the cell
counts as snow covered, the snow-free value otherwise. Implemented for every `AbstractSnowFraction`;
called from the `surface_exchange_coefficients!` kernel.
"""
function ground_roughness end

@inline function ground_roughness(::PointSnowFraction{Tf}, snow_roughness_scheme, i, j, state, surface) where {Tf}
    (; Ds) = state
    (; z0sf) = surface
    return column_sum(Ds, i, j) <= Tf(0.05) ? z0sf[i, j] : snow_roughness(snow_roughness_scheme, i, j, state, surface)
end

@inline function ground_roughness(::AbstractSnowFraction{Tf}, snow_roughness_scheme, i, j, state, surface) where {Tf}
    (; fsnow) = state
    (; z0sf) = surface
    return fsnow[i, j] <= eps(Tf) ? z0sf[i, j] : snow_roughness(snow_roughness_scheme, i, j, state, surface)
end

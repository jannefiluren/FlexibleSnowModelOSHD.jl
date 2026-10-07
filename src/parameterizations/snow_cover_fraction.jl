"""
$(TYPEDEF)

Binary snow cover fraction parameterization for a point location (snow or no snow). Unlike the
other schemes, ground roughness switches on a 5 cm snow-depth threshold and melt and sublimation
are not scaled up by `dfsnow_melt`.
"""
struct PointSnowFraction{Tf} <: AbstractSnowFraction{Tf} end

"""
$(TYPEDEF)

Snow cover fraction parameterization via a hyperbolic tangent of snow depth,
`fsnow = tanh(snowdepth / hfsn)`, adaptable via the single depth scale `hfsn`.

```jldoctest
using FlexibleSnowModelOSHD

TanhSnowFraction{Float32}(; hfsn = 0.3)

# output
TanhSnowFraction
├── hfsn: 0.3
└── dfsnow_melt: 0.25
```

# Fields

$(TYPEDFIELDS)
"""
@kwdef struct TanhSnowFraction{Tf} <: AbstractSnowFraction{Tf}
    "Snow-cover fraction depth scale (m)"
    hfsn::Tf = 0.1
    "Increase of fsnow for melt and sublimation, so thin snow does not deplete too slowly (-)"
    dfsnow_melt::Tf = 0.25
end

"""
$(TYPEDEF)

Seasonal snow cover fraction parameterization: the maximum of a seasonal SCF and a new-snow SCF.
The new-snow SCF is evaluated over two windows of the 14-day SWE history: the 14-day window
(since the 14-day SWE minimum) and the recent window (since the minimum preceding the most recent
SWE peak). Based on [Helbig et al. (2021)](https://doi.org/10.5194/tc-15-4607-2021) and
[Egli and Jonas (2009)](https://doi.org/10.1029/2008GL035545).

```jldoctest
using FlexibleSnowModelOSHD

SeasonalSnowFraction{Float32}(; nplateau = 3)

# output
SeasonalSnowFraction
├── sd_exp_hs: 0.549
├── sd_exp_slope: 0.309
├── sd_exp_flat: 0.84
├── c_season: 1.3
├── dswe_peak: 0.5
├── dswe_plateau: 0.1
├── nplateau: 3
├── a_nsnow: 0.2
├── b_nsnow: 0.3
├── dhs_min: 0.001
├── fsnow_min: 0.01
└── dfsnow_melt: 0.25
```

# Fields

$(TYPEDFIELDS)
"""
@kwdef struct SeasonalSnowFraction{Tf} <: AbstractSnowFraction{Tf}
    # Seasonal parameterization
    "Exponent of the seasonal maximum snow depth in the snow-depth standard deviation (Helbig et al., 2021) (-)"
    sd_exp_hs::Tf = 0.549
    "Exponent of the slope in the snow-depth standard deviation (Helbig et al., 2021) (-)"
    sd_exp_slope::Tf = 0.309
    "Exponent of the flat-field snow-depth standard deviation (Egli and Jonas, 2009), also used for new snow (-)"
    sd_exp_flat::Tf = 0.84
    "Coefficient of the seasonal SCF tanh (-)"
    c_season::Tf = 1.3

    # New snow parameterization (14-day and recent windows)
    "SWE drop ending the search for the most recent peak (kg/m^2)"
    dswe_peak::Tf = 0.5
    "SWE drop below which a day counts as plateau (kg/m^2)"
    dswe_plateau::Tf = 0.1
    "Consecutive plateau days ending the search for the preceding minimum"
    nplateau::Int = 2
    "Exponent on the maximum snow depth increase in the 14-day and recent new-snow SCF (-)"
    a_nsnow::Tf = 0.2
    "Scale of the 14-day and recent new-snow SCF (-)"
    b_nsnow::Tf = 0.3

    # Lower bounds
    "Minimum snow depth increase necessary to compute new snow SCF for numerical reasons (m)"
    dhs_min::Tf = 1.0e-3
    "Minimum SCF where snow is present (-)"
    fsnow_min::Tf = 0.01
    "Increase of fsnow for melt and sublimation, so thin snow does not deplete too slowly (-)"
    dfsnow_melt::Tf = 0.25
end

PointSnowFraction{Tf}(grid::Grid; kwargs...) where {Tf} = PointSnowFraction{Tf}()
SeasonalSnowFraction{Tf}(grid::Grid; kwargs...) where {Tf} = SeasonalSnowFraction{Tf}(; kwargs...)
TanhSnowFraction{Tf}(grid::Grid; kwargs...) where {Tf} = TanhSnowFraction{Tf}(; kwargs...)

"""
    ground_roughness(scheme, i, j, state, surface)

Roughness length of the ground at cell `(i, j)`: the snow value where the cell counts as
snow covered, the snow-free value otherwise. Implemented for every `AbstractSnowFraction`;
called from the `surface_exchange_coefficients!` kernel.
"""
function ground_roughness end

@inline function ground_roughness(::PointSnowFraction{Tf}, i, j, state, surface) where {Tf}
    (; Ds) = state
    (; z0_snow, z0sf) = surface
    return column_sum(Ds, i, j) <= Tf(0.05) ? z0sf[i, j] : z0_snow[i, j]
end

@inline function ground_roughness(::AbstractSnowFraction{Tf}, i, j, state, surface) where {Tf}
    (; fsnow) = state
    (; z0_snow, z0sf) = surface
    return fsnow[i, j] <= eps(Tf) ? z0sf[i, j] : z0_snow[i, j]
end

"""
    melt_snow_fraction(scheme, i, j, state)

Snow cover fraction used to scale melt and sublimation at cell `(i, j)`. Apart from the point
parameterization, `state.fsnow` is inflated by the scheme's `dfsnow_melt` to prevent that a
thin snow cover depletes unrealistically slowly. Any new scheme needs either a `dfsnow_melt`
field or its own method. Called from `snow!`.
"""
function melt_snow_fraction end

@inline function melt_snow_fraction(::PointSnowFraction{Tf}, i, j, state) where {Tf}
    (; fsnow) = state
    return fsnow[i, j]
end

@inline function melt_snow_fraction(scheme::AbstractSnowFraction{Tf}, i, j, state) where {Tf}
    (; fsnow) = state
    return min(fsnow[i, j] + scheme.dfsnow_melt, Tf(1.0))
end

"""
$(TYPEDSIGNATURES)

Snow cover fraction parameterizations for one grid cell.
"""
Base.@propagate_inbounds function snowcoverfraction_point!(
        scheme::AbstractSnowFraction, state, surface,
        snowdepth::Tf, SWEtmp::Tf, i::Integer, j::Integer, update_hist::Bool
    ) where {Tf <: Real}

    snow_covered_fraction!(scheme, state, surface, snowdepth, SWEtmp, i, j, update_hist)

    (; fsnow) = state
    # Final adjustments
    if snowdepth < eps(Tf)
        fsnow[i, j] = Tf(0.0)
    else
        fsnow[i, j] = min(fsnow[i, j], Tf(1.0))
    end

    return nothing
end

@inline function snow_covered_fraction!(
        ::PointSnowFraction{Tf}, state, surface,
        snowdepth::Tf, SWEtmp::Tf, i, j, update_hist::Bool
    ) where {Tf}
    (; fsnow) = state
    fsnow[i, j] = snowdepth > eps(Tf) ? Tf(1.0) : Tf(0.0)
    return nothing
end

@inline function snow_covered_fraction!(
        scheme::TanhSnowFraction{Tf}, state, surface,
        snowdepth::Tf, SWEtmp::Tf, i, j, update_hist::Bool
    ) where {Tf}
    (; fsnow) = state
    fsnow[i, j] = tanh(snowdepth / scheme.hfsn)
    return nothing
end

Base.@propagate_inbounds function snow_covered_fraction!(
        scheme::SeasonalSnowFraction{Tf}, state, surface,
        snowdepth::Tf, SWEtmp::Tf, i, j, update_hist::Bool
    ) where {Tf}
    (; sd_exp_hs, sd_exp_slope, sd_exp_flat, c_season) = scheme
    (; dswe_peak, dswe_plateau, nplateau, a_nsnow, b_nsnow, dhs_min, fsnow_min) = scheme
    (; fsnow, swehist, swemin, swemax, snowdepthhist, snowdepthmin, snowdepthmax) = state
    (; slopemu, xi, Ld) = surface

    # Today plus the 14-day SWE and snow-depth history (index 1 = today, 15 = oldest)
    SWEbuffer = MVector{15, Tf}(undef)
    snowdepthbuffer = MVector{15, Tf}(undef)
    SWEbuffer[1] = SWEtmp
    snowdepthbuffer[1] = snowdepth
    for k in 1:14
        SWEbuffer[k + 1] = swehist[k, i, j]
        snowdepthbuffer[k + 1] = snowdepthhist[k, i, j]
    end

    # 14-day window: SWE minimum of the buffer, and the maximum between today and that minimum
    imin_14d = first_argmin(SWEbuffer, 15)
    imax_14d = first_argmax(SWEbuffer, imin_14d)

    # Recent window: most recent SWE peak, walking back while SWE drops by less than dswe_peak
    imax_recent = 1
    for k in 1:14
        if SWEbuffer[k + 1] - SWEbuffer[k] >= -dswe_peak
            imax_recent = k + 1
        else
            break
        end
    end

    # ... and the minimum preceding that peak; nplateau days with drop < dswe_plateau end the search
    imin_recent = imax_recent
    if imax_recent < 15
        imin_recent = imax_recent + 1
        plateau_days = 0
        for k in (imax_recent + 1):14
            if SWEbuffer[k + 1] < SWEbuffer[imin_recent]
                imin_recent = k + 1
            end
            plateau_days = SWEbuffer[k + 1] - SWEbuffer[k] > -dswe_plateau ? plateau_days + 1 : 0
            if plateau_days >= nplateau
                break
            end
        end
    end

    # manual loop: a view-based maximum on the MVector scratch risks heap allocation
    snowdepthmax_recent = snowdepthbuffer[imax_recent]
    for k in (imax_recent + 1):imin_recent
        snowdepthmax_recent = max(snowdepthmax_recent, snowdepthbuffer[k])
    end

    # New snow stored on old snow: current and maximum snow depth increase in both windows
    dsnowdepth_14d = snowdepth - snowdepthbuffer[imin_14d]
    if (dsnowdepth_14d < eps(Tf))
        dsnowdepth_14d = Tf(0)
    end
    dsnowdepth_14d_max = snowdepthbuffer[imax_14d] - snowdepthbuffer[imin_14d]
    if (dsnowdepth_14d_max < eps(Tf))
        dsnowdepth_14d_max = Tf(0)
    end

    dsnowdepth_recent = snowdepth - snowdepthbuffer[imin_recent]
    if (dsnowdepth_recent < eps(Tf))
        dsnowdepth_recent = Tf(0)
    end
    dsnowdepth_recent_max = snowdepthmax_recent - snowdepthbuffer[imin_recent]
    if (dsnowdepth_recent_max < eps(Tf))
        dsnowdepth_recent_max = Tf(0)
    end

    # a max snow depth increase below the current increase would inflate fsnow
    dsnowdepth_14d_max = max(dsnowdepth_14d_max, dsnowdepth_14d)
    dsnowdepth_recent_max = max(dsnowdepth_recent_max, dsnowdepth_recent)

    # Season-long SWE and snow-depth extremes; reset when snow-free, min restarts at each new max
    if (SWEtmp < eps(Tf))
        swemax[i, j] = Tf(0)
        swemin[i, j] = Tf(0)
    end
    if (snowdepth < eps(Tf))
        snowdepthmax[i, j] = Tf(0)
        snowdepthmin[i, j] = Tf(0)
    end

    if (SWEtmp >= swemax[i, j])
        swemax[i, j] = SWEtmp
        swemin[i, j] = SWEtmp
    end
    # snowdepth can exceed snowdepthmax since the max index is chosen on SWE
    if (snowdepth >= snowdepthmax[i, j])
        snowdepthmax[i, j] = snowdepth
        snowdepthmin[i, j] = snowdepth
    end

    if (SWEtmp < swemax[i, j] && SWEtmp < swemin[i, j])
        swemin[i, j] = SWEtmp
    end
    if (snowdepth < snowdepthmax[i, j] && snowdepth < snowdepthmin[i, j])
        snowdepthmin[i, j] = snowdepth
    end

    # Seasonal SCF (Helbig et al., 2021) with topography-dependent snow-depth standard deviation
    fsnow_season = Tf(0)
    sd_snowdepth1 = exp(Tf(-1) / (Ld[i, j] / xi[i, j])^Tf(2))
    sd_snowdepth2 = snowdepthmax[i, j]^sd_exp_hs
    sd_snowdepth3 = slopemu[i, j]^sd_exp_slope
    sd_snowdepth0 = sd_snowdepth1 * sd_snowdepth2 * sd_snowdepth3
    # Flat pixels use a slope-free standard deviation (Egli and Jonas, 2009)
    if (!(slopemu[i, j] > eps(Tf)))
        sd_snowdepth0 = snowdepthmax[i, j]^sd_exp_flat
    end
    if (snowdepthmax[i, j] > eps(Tf))
        fsnow_season = tanh(c_season * snowdepthmin[i, j] / sd_snowdepth0)
    end

    # New-snow SCF in the 14-day and recent windows (flat-field standard deviation)
    fsnow_new_14d = Tf(0)
    if (dsnowdepth_14d > eps(Tf) && dsnowdepth_14d_max > dhs_min)
        sd_snowdepth0_14d = dsnowdepth_14d_max^sd_exp_flat
        fsnow_new_14d = tanh(dsnowdepth_14d / sd_snowdepth0_14d + dsnowdepth_14d / dsnowdepth_14d_max^a_nsnow / b_nsnow)
    end

    fsnow_new_recent = Tf(0)
    if (dsnowdepth_recent > eps(Tf) && dsnowdepth_recent_max > dhs_min)
        sd_snowdepth0_recent = dsnowdepth_recent_max^sd_exp_flat
        fsnow_new_recent = tanh(dsnowdepth_recent / sd_snowdepth0_recent + dsnowdepth_recent / dsnowdepth_recent_max^a_nsnow / b_nsnow)
    end

    fsnow[i, j] = max(fsnow_season, fsnow_new_14d, fsnow_new_recent, fsnow_min)

    # Roll the 14-day SWE/depth history
    if update_hist
        for k in 1:14
            swehist[k, i, j] = SWEbuffer[k]
            snowdepthhist[k, i, j] = snowdepthbuffer[k]
        end
    end

    return nothing
end

"""
$(TYPEDSIGNATURES)

Snow cover fraction for one grid cell (host convenience wrapper around
[`snowcoverfraction_point!`](@ref), kept for API compatibility).

The buffer arguments are accepted but ignored: the history buffers are now
function-local (they were always pure workspace).
"""
function snow_cover_fraction!(fsm::FSM{Tf}, snowdepth::Tf, SWEtmp::Tf, t::DateTime, i::Int, j::Int, SWEbuffer::AbstractArray{Tf}, snowdepthbuffer::AbstractArray{Tf}, diffSWEbuffer::AbstractArray{Tf}) where {Tf <: Real}

    # update history of SWE and hs only if they correspond to 6:00am values
    update_hist = 4.5 < hour(t) < 5.5

    @inbounds snowcoverfraction_point!(
        fsm.physics.snow_fraction, fsm.state, fsm.surface,
        snowdepth, SWEtmp, i, j, update_hist
    )

    return nothing
end

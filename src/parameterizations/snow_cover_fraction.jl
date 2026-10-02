struct SeasonalSnowFraction{Tf} <: AbstractSnowFraction{Tf} end
struct HelbigSnowFraction{Tf} <: AbstractSnowFraction{Tf} end
struct HelbigMaxSnowFraction{Tf} <: AbstractSnowFraction{Tf} end
struct PointSnowFraction{Tf} <: AbstractSnowFraction{Tf} end
@kwdef struct TanhSnowFraction{Tf} <: AbstractSnowFraction{Tf}
    hfsn::Tf = 0.1             # Snow-cover fraction depth scale (m)
end

SeasonalSnowFraction{Tf}(grid::Grid; kwargs...) where {Tf} = SeasonalSnowFraction{Tf}()
HelbigSnowFraction{Tf}(grid::Grid; kwargs...) where {Tf} = HelbigSnowFraction{Tf}()
HelbigMaxSnowFraction{Tf}(grid::Grid; kwargs...) where {Tf} = HelbigMaxSnowFraction{Tf}()
PointSnowFraction{Tf}(grid::Grid; kwargs...) where {Tf} = PointSnowFraction{Tf}()
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
model, `state.fsnow` is inflated to prevent that a thin snow cover deplets unrealistically
slow. Implemented for every `AbstractSnowFraction`; called from the `snow!`.
"""
function melt_snow_fraction end

@inline function melt_snow_fraction(::PointSnowFraction{Tf}, i, j, state) where {Tf}
    (; fsnow) = state
    return fsnow[i, j]
end

@inline function melt_snow_fraction(::AbstractSnowFraction{Tf}, i, j, state) where {Tf}
    (; fsnow) = state
    return min(fsnow[i, j] + Tf(0.25), Tf(1.0))
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

Base.@propagate_inbounds function snow_covered_fraction!(
        ::SeasonalSnowFraction{Tf}, state, surface,
        snowdepth::Tf, SWEtmp::Tf, i, j, update_hist::Bool
    ) where {Tf}
    (; fsnow, swehist, swemin, swemax, snowdepthhist, snowdepthmin, snowdepthmax) = state
    (; slopemu, xi, Ld) = surface
    # topo terms for the snow-depth standard deviation
    sd_snowdepth1 = exp(Tf(-1) / (Ld[i, j] / xi[i, j])^Tf(2))
    sd_snowdepth3 = slopemu[i, j]^Tf(0.309)

    # merge current SWEtmp with SWEtmp history from past 14 days
    SWEbuffer = MVector{15, Tf}(undef)
    snowdepthbuffer = MVector{15, Tf}(undef)
    SWEbuffer[1] = SWEtmp
    snowdepthbuffer[1] = snowdepth
    for k in 1:14
        SWEbuffer[k + 1] = swehist[k, i, j]
        snowdepthbuffer[k + 1] = snowdepthhist[k, i, j]
    end

    # indices of the global min and max in the SWE buffer
    iabsmin = first_argmin(SWEbuffer, 15)

    # maximum between today and the minimum (inclusive)
    # iabsmax = first_argmax(@view SWEbuffer[1:iabsmin], iabsmin)
    iabsmax = first_argmax(view(SWEbuffer, 1:iabsmin), iabsmin)


    # -------------------------------------------------------------
    # Find recent SWE maximum (today -> older)
    # -------------------------------------------------------------
    maxidx = 1

    for k in 1:14
        diffSWEbuffer = SWEbuffer[k + 1] - SWEbuffer[k]

        # Continue while SWE is increasing or nearly flat
        if diffSWEbuffer >= Tf(-0.5)
            maxidx = k + 1
        else
            break
        end
    end

    # -------------------------------------------------------------
    # Find preceding SWE minimum with plateau hysteresis
    # -------------------------------------------------------------
    if maxidx == 15

        irecentmin = 15
        recent_max = snowdepthbuffer[15]

    else

        irecentmin = maxidx
        plateau_days = 0
        started = false

        for k in maxidx:14

            diffSWEbuffer = SWEbuffer[k + 1] - SWEbuffer[k]

            if !started

                # First meaningful decrease starts the event
                if diffSWEbuffer < Tf(-0.5)  # not considering 0.5mm SWE gain a new snowfall TZ 10.26
                    started = true

                    if SWEbuffer[k + 1] < SWEbuffer[irecentmin]
                        irecentmin = k + 1
                    end
                end

            else

                # Update running minimum
                if SWEbuffer[k + 1] < SWEbuffer[irecentmin]
                    irecentmin = k + 1
                end

                # Plateau detection
                if diffSWEbuffer > Tf(-0.1)
                    plateau_days += 1
                else
                    plateau_days = 0
                end

                if plateau_days >= 2    #one day of no increase does not stop the recent calculation TZ 10.26
                    break
                end
            end
        end

        recent_max = maximum(@view snowdepthbuffer[maxidx:irecentmin])

    end

    # use indices to determine snowdepth amounts
    snowdepthmin_buffer = snowdepthbuffer[iabsmin]
    snowdepthmax_buffer = snowdepthbuffer[iabsmax]
    snowdepthmin_recent = snowdepthbuffer[irecentmin]
    snowdepthmax_recent = recent_max

    # Compute storage of new snow on old snow in snowdepthbuffer
    dsnowdepth = snowdepth - snowdepthmin_buffer
    if (dsnowdepth < eps(Tf))
        dsnowdepth = Tf(0)
    end

    # compute maximum new snow depth (dswemax) in SWEbuffer
    dsnowdepthmax = snowdepthmax_buffer - snowdepthmin_buffer
    if (dsnowdepthmax < eps(Tf))
        dsnowdepthmax = Tf(0)
    end

    # cap dsnowdepthmax at dsnowdepth, else fsnow inflates
    if (dsnowdepthmax < dsnowdepth)
        dsnowdepthmax = dsnowdepth
    end

    # recent new snow stored on old snow
    dsnowdepth_recent = snowdepth - snowdepthmin_recent
    if (dsnowdepth_recent < eps(Tf))
        dsnowdepth_recent = Tf(0)
    end

    dsnowdepth_recent_max = max.(snowdepthmax_recent .- snowdepthmin_recent, 0)
    if (dsnowdepth_recent_max < eps(Tf))
        dsnowdepth_recent_max = Tf(0)
    end
    # don't accept dsnowdepthmax to be smaller then dsnowdepth, otherwise larger fnsnow values
    if (dsnowdepth_recent_max < dsnowdepth_recent)
        dsnowdepth_recent_max = dsnowdepth_recent
    end
    # state variables interpreting the whole SWEtmp history, not only the past 14 days in the buffer
    # Set swemax and swemin equal to zero if no snow, same with corresponding snow depth values
    if (SWEtmp < eps(Tf))
        swemax[i, j] = Tf(0)
        swemin[i, j] = Tf(0)
    end
    if (snowdepth < eps(Tf))
        snowdepthmax[i, j] = Tf(0)
        snowdepthmin[i, j] = Tf(0)
    end

    # Set swemax and swemin equal to SWEtmp if maximum, store also snowdepthmax and snowdepthmin of those time steps
    if (SWEtmp >= swemax[i, j])
        swemax[i, j] = SWEtmp
        swemin[i, j] = SWEtmp
    end

    # snowdepth can exceed snowdepthmax since the max index is chosen on SWE
    if (snowdepth >= snowdepthmax[i, j])
        snowdepthmax[i, j] = snowdepth
        snowdepthmin[i, j] = snowdepth
    end

    # Set swemin equal SWEtmp if smaller than swemin, same with corresponding snow depth value
    if (SWEtmp < swemax[i, j] && SWEtmp < swemin[i, j])
        swemin[i, j] = SWEtmp
    end
    if (snowdepth < snowdepthmax[i, j] && snowdepth < snowdepthmin[i, j])
        snowdepthmin[i, j] = snowdepth
    end

    # Initial guess of snow covered fraction
    fsnow_season = Tf(0)

    # Seasonal SCF (Helbig et al.; Egli & Jonas)
    # standard deviation
    sd_snowdepth2 = snowdepthmax[i, j]^Tf(0.549)
    sd_snowdepth0 = sd_snowdepth1 * sd_snowdepth2 * sd_snowdepth3
    # flat pixels use a slope-free standard deviation (Egli)
    if (!(slopemu[i, j] > eps(Tf)))
        sd_snowdepth0 = snowdepthmax[i, j]^Tf(0.84)
    end
    if (snowdepthmax[i, j] > eps(Tf))
        fsnow_season = tanh(Tf(1.3) * snowdepthmin[i, j] / sd_snowdepth0)
    end

    coeff_vari = sd_snowdepth0 / snowdepthmax[i, j]

    # SCF from the last 14 days' new snow (flat-field standard deviation)
    fsnow_nsnow = Tf(0)

    sd_snowdepth0_dhs = dsnowdepthmax^Tf(0.84)
    if (dsnowdepthmax > eps(Tf))
        fsnow_nsnow = tanh(dsnowdepth / sd_snowdepth0_dhs + dsnowdepth / dsnowdepthmax^0.2 / Tf(0.3))
    end

    # SCF from new snow since the last minimum (flat-field standard deviation)
    fsnow_nsnow_recent = Tf(0)

    sd_snowdepth0_dhs_recent = dsnowdepth_recent^Tf(0.84)
    if (dsnowdepth_recent > eps(Tf)&& (dsnowdepth_recent_max > Tf(1.0e-3)))  #meaningful snowheight due to swee search
        fsnow_nsnow_recent = tanh(dsnowdepth_recent / sd_snowdepth0_dhs_recent + dsnowdepth_recent / dsnowdepth_recent_max^0.2 / Tf(0.3))
    end

    fsnow_nsnow = max(fsnow_nsnow, fsnow_nsnow_recent)


    fsnow[i, j] = max(fsnow_season, fsnow_nsnow)

    # refresh the 14-day SWE/depth history
    if update_hist
        for k in 1:14
            swehist[k, i, j] = SWEbuffer[k]
            snowdepthhist[k, i, j] = snowdepthbuffer[k]
        end
    end

    fsnow[i, j] = max(fsnow[i, j], Tf(0.01))
    return nothing
end

@inline function snow_covered_fraction!(
        ::HelbigSnowFraction{Tf}, state, surface,
        snowdepth::Tf, SWEtmp::Tf, i, j, update_hist::Bool
    ) where {Tf}
    (; fsnow) = state
    (; slopemu, xi, Ld) = surface
    sd_snowdepth2 = snowdepth^Tf(0.549)
    sd_snowdepth1 = exp(Tf(-1) / (Ld[i, j] / xi[i, j])^Tf(2))
    sd_snowdepth3 = slopemu[i, j]^Tf(0.309)
    sd_snowdepth0 = sd_snowdepth1 * sd_snowdepth2 * sd_snowdepth3

    fsnow[i, j] = tanh(Tf(1.3) * snowdepth / sd_snowdepth0)
    return nothing
end

@inline function snow_covered_fraction!(
        ::HelbigMaxSnowFraction{Tf}, state, surface,
        snowdepth::Tf, SWEtmp::Tf, i, j, update_hist::Bool
    ) where {Tf}
    (; fsnow, snowdepthmax) = state
    (; slopemu, xi, Ld) = surface
    if snowdepth == Tf(0)
        snowdepthmax[i, j] = Tf(0.0)
    end

    if snowdepth > snowdepthmax[i, j]
        snowdepthmax[i, j] = snowdepth
    end

    sd_snowdepth2 = snowdepthmax[i, j]^Tf(0.549)
    sd_snowdepth1 = exp(Tf(-1) / (Ld[i, j] / xi[i, j])^Tf(2))
    sd_snowdepth3 = slopemu[i, j]^Tf(0.309)
    sd_snowdepth0 = sd_snowdepth1 * sd_snowdepth2 * sd_snowdepth3

    fsnow[i, j] = tanh(Tf(1.3) * snowdepth / sd_snowdepth0)
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

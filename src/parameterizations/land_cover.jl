# Land-cover scheme: OpenCover (open terrain / glacier) and ForestCover (one-layer canopy), together
# with all AbstractLandCover-dispatched physics — canopy mass balance, surface-layer turbulent
# exchange, shortwave/longwave radiative transfer, and the surface (and canopy) energy balance. The
# process files (canopy.jl, surface_exchange_coefficients.jl, radiation.jl, surface_energy_balance.jl)
# are thin launchers that call these methods.

struct OpenCover{Tf, S <: AbstractStabilityCorrection{Tf}} <: AbstractLandCover{Tf}
    # Atmospheric stability correction for the open-terrain eddy diffusivity. A type parameter
    # (so the scheme stays isbits); its default is on the constructor, since S follows the value.
    stability::S
end

@kwdef struct ForestCover{Tf} <: AbstractLandCover{Tf}
    fsar::Tf = 0.1                       # Albedo adjustment range vs vegetation fraction (-)
    avg0::Tf = 0.1                       # Snow-free vegetation albedo (-)
    avgs::Tf = 0.4                       # Snow-covered vegetation albedo (-)
    psf::Tf = 1                          # Solid precipitation multiplier at min canopy cover (-)
    psr::Tf = 0.1                        # Solid precipitation multiplier range (-)
    tcnc::Tf = 3600 * 240                # Canopy unloading time scale for cold snow (s)
    tcnm::Tf = 3600 * 48                 # Canopy unloading time scale for melting snow (s)
    rchd::Tf = 0.67                      # Ratio of displacement height to canopy height (-)
    rchz::Tf = 0.2                       # Ratio of roughness length to canopy height (-)
    zgf::Tf = 1                          # Roughness length adjustment factor vs vegetation fraction (-)
    zgr::Tf = 0                          # Roughness length adjustment range vs vegetation fraction (-)
    wcan::Tf = 2.5                       # Parameter of exponential wind profile (-)
    khcf::Tf = 3                         # Diffusivity adjustment for canopy effects (-)
    cveg::Tf = 20                        # Vegetation turbulent transfer coefficient ((s/m)^0.5)
    gsnf::Tf = 0                         # Snow-free vegetation moisture conductance (m/s)
    zsub::Tf = 2                         # Sub-canopy reference height (m)
end

# Default to the Louis (1982) stability correction over open terrain, the realistic choice used
# across the model; pass stability = NoStabilityCorrection{Tf}() for a neutral surface layer.
OpenCover{Tf}(; stability = LouisStabilityCorrection{Tf}()) where {Tf} =
    OpenCover{Tf, typeof(stability)}(stability)
OpenCover{Tf}(grid::Grid; kwargs...) where {Tf} = OpenCover{Tf}(; kwargs...)
ForestCover{Tf}(grid::Grid; kwargs...) where {Tf} = ForestCover{Tf}(; kwargs...)

# Per-cell mix of open and forest land cover for a single run. `forestcells` is a per-cell mask:
# physics uses the forest scheme where it is true and the open scheme elsewhere (glacier cells use
# open cover). Carries a grid-sized array like the per-cell albedo schemes, so it is `@adapt_structure`d.
struct MixedLandCover{Tf, O, F, MB} <: AbstractLandCover{Tf}
    open::O
    forest::F
    forestcells::MB
end
MixedLandCover(open::OpenCover{Tf}, forest::ForestCover{Tf}, forestcells::AbstractMatrix{Bool}) where {Tf} =
    MixedLandCover{Tf, typeof(open), typeof(forest), typeof(forestcells)}(open, forest, forestcells)
@adapt_structure MixedLandCover

canopy_fsar(c::ForestCover) = c.fsar
canopy_avg0(c::ForestCover) = c.avg0
canopy_avgs(c::ForestCover) = c.avgs

"""
    canopy_snow!(land_cover, i, j, state, diag, surface, params)

Snow on the canopy at cell `(i, j)`: interception from the throughfall `diag.Sfeff`,
sublimation and unloading, updating `state.Sveg` and `diag.intcpt`/`Sbveg`/`unload`.
"""
function canopy_snow! end

@inline canopy_snow!(::OpenCover, i, j, state, diag, surface, params) = nothing

@inline function canopy_snow!(land_cover::ForestCover{Tf}, i, j, state, diag, surface, params) where {Tf}

    @unpack_constants(Tf)

    (; dt) = params
    (; tcnc, tcnm) = land_cover
    (; scap, fveg, pmultf) = surface
    (; Sveg, Tveg) = state
    (; unload, intcpt, Sbveg, Sfeff, Eveg) = diag

    unload[i, j] = Tf(0)
    intcpt[i, j] = Tf(0)
    Sbveg[i, j] = Tf(0)

    # Remove precipitation scaling applied to forcing data
    Sfeff[i, j] = pmultf[i, j] * Sfeff[i, j]

    # Interception of snow on canopies
    intcpt[i, j] = (scap[i, j] - Sveg[i, j]) * (Tf(1) - exp(-fveg[i, j] * Sfeff[i, j] * dt / scap[i, j]))
    Sveg[i, j] = Sveg[i, j] + intcpt[i, j]
    Sfeff[i, j] = Sfeff[i, j] - intcpt[i, j] / dt

    # Preferential deposition of snowfall in canopy gaps (not mass conserving)
    Sfeff[i, j] = (land_cover.psf - land_cover.psr * fveg[i, j]) * Sfeff[i, j]

    # Sublimation of intercepted snow
    Evegs = Tf(0)
    if (Sveg[i, j] > eps(Tf) || Tveg[i, j] < Tm)
        Evegs = Eveg[i, j]
    end
    Sveg[i, j] = Sveg[i, j] - Evegs * dt
    Sbveg[i, j] = Evegs * dt
    if (Sveg[i, j] < Tf(0))
        Sbveg[i, j] = Sbveg[i, j] + Sveg[i, j]
    end
    Sveg[i, j] = max(Sveg[i, j], Tf(0))

    # Unloading of intercepted snow
    tunl = tcnc
    if (Tveg[i, j] >= Tm)
        tunl = tcnm
    end
    tunl = max(tunl, dt)
    unload[i, j] = Sveg[i, j] * dt / tunl
    Sveg[i, j] = Sveg[i, j] - unload[i, j]

    return nothing
end


"""
    exchange_coefficients!(land_cover, i, j, state, diag, surface, params, meteo, z0g)

Eddy diffusivities for turbulent heat and moisture transfer at cell `(i, j)`, implemented
for every `AbstractLandCover`: `diag.KH` and `diag.KWg` over open terrain, and
`diag.KHa`, `diag.KHg`, `diag.KHv`, `diag.KWg`, `diag.KWv` under a canopy. The ground
roughness length `z0g` is resolved by the caller; the wind and temperature reference
heights are derived from the measurement heights in `params` inside each method — open
terrain uses them directly, a canopy offsets them by the canopy height.
"""
function exchange_coefficients! end

@inline function exchange_coefficients!(c::OpenCover{Tf}, i, j, state, diag, surface, params, meteo, z0g) where {Tf}
    @unpack_constants(Tf)
    (; zU, zT) = params
    (; Sice, Tsrf) = state
    (; KH, KWg, gs1, Qa, Uaeff) = diag
    (; Ta, Ps) = meteo

    # Roughness lengths and friction velocity
    z0 = z0g
    z0h = Tf(0.1) * z0
    CD = (vkman / log(zU / z0))^Tf(2)
    ustar = sqrt(CD) * Uaeff[i, j]

    fh = stability_factor(c.stability, CD, z0, Ta[i, j], Tsrf[i, j], Uaeff[i, j], zU, zT)

    # Eddy diffusivities
    KH[i, j] = fh * vkman * ustar / log(zT / z0h)
    Qs = qsat(Ps[i, j], Tsrf[i, j])
    if (Sice[1, i, j] > eps(Tf) || Qa[i, j] > Qs)
        KWg[i, j] = KH[i, j]
    else
        KWg[i, j] = gs1[i, j] * KH[i, j] / (gs1[i, j] + KH[i, j])
    end
    return nothing
end

@inline function exchange_coefficients!(c::ForestCover{Tf}, i, j, state, diag, surface, params, meteo, z0g) where {Tf}
    @unpack_constants(Tf)
    (; zU, zT) = params
    (; zsub, gsnf) = c
    (; fveg, fves, VAI, hcan) = surface
    (; Sveg, Tsrf, Tveg, Qcan) = state
    (; KHa, KHg, KHv, KWg, KWv, Usc, gs1, Uaeff) = diag
    (; Ps) = meteo

    # Reference heights measured above the canopy
    zU1 = zU + hcan[i, j]
    zT1 = zT + hcan[i, j]

    # Roughness lengths, friction velocity and canopy wind profile
    z0g = (c.zgf + c.zgr * fveg[i, j]) * z0g
    z0h = Tf(0.1) * z0g
    dh = c.rchd * hcan[i, j]
    z0v = c.rchz * hcan[i, j]
    ustar = vkman * Uaeff[i, j] / log((zU1 - dh) / z0v)
    Uh = (ustar / vkman) * log((hcan[i, j] - dh) / z0v)
    KHh = vkman * ustar * (hcan[i, j] - dh)
    Usf = exp(c.wcan * (zsub / hcan[i, j] - Tf(1))) * Uh

    Uso = Uaeff[i, j] * log(zsub / z0g) / log(zU / z0g)

    # Eddy diffusivities
    rad = (log((zT1 - dh) / (hcan[i, j] - dh)) / (vkman * ustar) + hcan[i, j] * (exp(c.wcan * (Tf(1) - (z0v + dh) / hcan[i, j])) - Tf(1)) / (c.wcan * KHh)) / c.khcf
    KHa[i, j] = sqrt(fves[i, j]) / rad
    Usub = sqrt(fves[i, j]) * Usf + (Tf(1) - sqrt(fves[i, j])) * Uso
    Usub = max(Usub, Tf(0.1))
    rgd = Tf(1) / (vkman^Tf(2) * Usub) * log(zsub / z0h) * log(zsub / z0g)
    KHg[i, j] = Tf(1) / rgd
    Uc = exp(c.wcan * ((z0v + dh) / hcan[i, j] - Tf(1))) * Uh
    KHv[i, j] = VAI[i, j] * sqrt(Uc) / c.cveg
    # Usc leaves the model only through the OSHDinternal output catalog (uaca)
    Usc[i, j] = Usub

    Qs = qsat(Ps[i, j], Tsrf[i, j])
    if (Qcan[i, j] > Qs)
        KWg[i, j] = KHg[i, j]
    else
        KWg[i, j] = gs1[i, j] * KHg[i, j] / (gs1[i, j] + KHg[i, j])
    end
    Qs = qsat(Ps[i, j], Tveg[i, j])
    if (Sveg[i, j] > eps(Tf) || Qcan[i, j] > Qs)
        KWv[i, j] = KHv[i, j]
    else
        KWv[i, j] = gsnf * KHv[i, j] / (gsnf + KHv[i, j])
    end
    return nothing
end


"""
    solar_radiation!(land_cover, i, j, state, diag, surface, meteo)

Surface albedo and shortwave transmission for cell `(i, j)`: fills
`diag.SWsrf`, `diag.SWveg` and `diag.SWsci`. Expects `state.albs`
to already hold the bare-ground albedo where the snow has gone.
"""
function solar_radiation! end

@inline function solar_radiation!(c::OpenCover{Tf}, i, j, state, diag, surface, meteo) where {Tf}
    (; albs) = state
    (; SWveg, SWsrf, SWsci) = diag
    (; Sdif, Sdir) = meteo

    asrf = albs[i, j]
    SWveg[i, j] = Tf(0)
    SWsrf[i, j] = (Tf(1) - asrf) * (Sdir[i, j] + Sdif[i, j])
    SWsci[i, j] = Sdif[i, j] + Sdir[i, j]
    return nothing
end

@inline function solar_radiation!(c::ForestCover{Tf}, i, j, state, diag, surface, meteo) where {Tf}
    (; albs, fsnow, Sveg) = state
    (; SWveg, SWsrf, SWsci) = diag
    (; fveg, fsky, fsky_terr, scap, trcn) = surface
    (; Sdif, Sdir, Tv) = meteo

    asrf = albs[i, j]
    if (fsnow[i, j] > eps(Tf))
        asrf *= Tf(1) - fveg[i, j] * canopy_fsar(c)
    end

    fcans = Tf(0.0)
    if (scap[i, j] > eps(Tf))
        fcans = Sveg[i, j] / scap[i, j]
    end
    aveg = (Tf(1) - fcans) * canopy_avg0(c) + fcans * canopy_avgs(c)

    Sdif_aux = fsky[i, j] / fsky_terr[i, j] * Sdif[i, j]
    tdif = trcn[i, j]
    tdir = Tv[i, j]
    SWsrf[i, j] = (Tf(1) - asrf) * (tdif * Sdif_aux + tdir * Sdir[i, j])
    SWveg[i, j] = ((Tf(1) - tdif) * (Tf(1) - aveg) + tdif * asrf * (Tf(1) - tdif)) * Sdif_aux + (tdir * fveg[i, j] * (Tf(1) - aveg) + tdir * asrf * (Tf(1) - tdif)) * Sdir[i, j]
    SWsci[i, j] = tdif * Sdif_aux + tdir * Sdir[i, j]
    return nothing
end

"""
    thermal_radiation!(land_cover, i, j, diag, surface, meteo)

Effective incoming longwave `diag.LWeff` for cell `(i, j)`. Without canopy the terrain
emission is computed here, while in forested cells `energy_balance!` accounts for it.
"""
function thermal_radiation! end

@inline function thermal_radiation!(c::OpenCover{Tf}, i, j, diag, surface, meteo) where {Tf}
    @unpack_constants(Tf)
    (; LWeff) = diag
    (; fsky_terr) = surface
    (; LW, Ta) = meteo

    LWeff[i, j] = fsky_terr[i, j] * LW[i, j] + (Tf(1) - fsky_terr[i, j]) * sb * Ta[i, j]^Tf(4)
    return nothing
end

@inline function thermal_radiation!(c::ForestCover{Tf}, i, j, diag, surface, meteo) where {Tf}
    @unpack_constants(Tf)
    (; LWeff) = diag
    (; fsky) = surface
    (; Ta, LW) = meteo

    LWeff[i, j] = fsky[i, j] * LW[i, j] + (Tf(1) - fsky[i, j]) * sb * Ta[i, j]^Tf(4)
    return nothing
end

"""
    energy_balance!(land_cover, substrate, i, j, state, diag, surface, params, meteo)

Solve the energy balance at cell `(i, j)`, implemented for every `AbstractLandCover`:
`OpenCover` solves the surface alone, `ForestCover` solves the joint surface and canopy
system. `substrate` caps the surface temperature over bare glacier ice and is
unused under a canopy.
"""
function energy_balance! end

@inline function energy_balance!(::OpenCover{Tf}, substrate, i, j, state, diag, surface, params, meteo) where {Tf}

    @unpack_constants(Tf)

    (; dt) = params
    (; trcn) = surface
    (; Sice, Tcan, Tsrf, Tveg) = state
    (; Esrf, G, H, Hsrf, LE, LEsrf, LWsci, LWveg, Melt, Rnet, Rsrf, SWsrf, Ds1, Ts1, ks1, KH, KWg, Qa, LWeff) = diag
    (; Ps, Ta) = meteo

    # Reported as air temperature so that Tcan and Tveg are defined in open runs
    Tveg[i, j] = Ta[i, j]
    Tcan[i, j] = Ta[i, j]

    # Saturation humidity and density of air
    Qs = qsat(Ps[i, j], Tsrf[i, j])
    Lh = Lv
    if (Tsrf[i, j] < Tm || Sice[1, i, j] > eps(Tf))
        Lh = Ls
    end
    D = Lh * Qs / (Rwat * Tsrf[i, j]^Tf(2))
    rho = Ps[i, j] / (Rair * Ta[i, j])

    # Explicit fluxes
    Esrf[i, j] = rho * KWg[i, j] * (Qs - Qa[i, j])
    G[i, j] = Tf(2) * ks1[i, j] * (Tsrf[i, j] - Ts1[i, j]) / Ds1[i, j]
    H[i, j] = cp * rho * KH[i, j] * (Tsrf[i, j] - Ta[i, j])
    LE[i, j] = Lh * Esrf[i, j]
    Melt[i, j] = Tf(0)
    Rnet[i, j] = SWsrf[i, j] + trcn[i, j] * LWeff[i, j] - sb * Tsrf[i, j]^Tf(4) + (Tf(1) - trcn[i, j]) * sb * Tveg[i, j]^Tf(4)

    # Surface energy balance increments without melt
    dTs = (Rnet[i, j] - G[i, j] - H[i, j] - LE[i, j]) / (Tf(4) * sb * Tsrf[i, j]^Tf(3) + Tf(2) * ks1[i, j] / Ds1[i, j] + rho * (cp * KH[i, j] + Lh * D * KWg[i, j]))
    dE = rho * KWg[i, j] * D * dTs
    dG = Tf(2) * ks1[i, j] * dTs / Ds1[i, j]
    dH = cp * rho * KH[i, j] * dTs
    dR = Tf(-4) * sb * Tsrf[i, j]^Tf(3) * dTs

    # Surface melting
    if (Tsrf[i, j] + dTs > Tm && Sice[1, i, j] > eps(Tf))
        Melt[i, j] = column_sum(Sice, i, j) / dt
        dTs = (Rnet[i, j] - G[i, j] - H[i, j] - LE[i, j] - Lf * Melt[i, j]) / (Tf(4) * sb * Tsrf[i, j]^Tf(3) + Tf(2) * ks1[i, j] / Ds1[i, j] + rho * (cp * KH[i, j] + Ls * D * KWg[i, j]))
        dE = rho * KWg[i, j] * D * dTs
        dG = Tf(2) * ks1[i, j] * dTs / Ds1[i, j]
        dH = cp * rho * KH[i, j] * dTs
        dR = Tf(-4) * sb * Tsrf[i, j]^Tf(3) * dTs
        if (Tsrf[i, j] + dTs < Tm)
            Qs = qsat(Ps[i, j], Tm)
            Esrf[i, j] = rho * KWg[i, j] * (Qs - Qa[i, j])
            G[i, j] = Tf(2) * ks1[i, j] * (Tm - Ts1[i, j]) / Ds1[i, j]
            H[i, j] = cp * rho * KH[i, j] * (Tm - Ta[i, j])
            LE[i, j] = Ls * Esrf[i, j]
            Rnet[i, j] = SWsrf[i, j] + trcn[i, j] * LWeff[i, j] - sb * Tm^Tf(4) + (Tf(1) - trcn[i, j]) * sb * Tveg[i, j]^Tf(4)
            Melt[i, j] = (Rnet[i, j] - H[i, j] - LE[i, j] - G[i, j]) / Lf
            Melt[i, j] = max(Melt[i, j], Tf(0.0))
            dE = Tf(0.0)
            dG = Tf(0.0)
            dH = Tf(0.0)
            dR = Tf(0.0)
            dTs = Tm - Tsrf[i, j]
        end
    end

    # Bare glacier ice as an infinite heat reservoir: capping Tsrf at melting point is not energy-conserving
    if substrate isa IceSubstrate
        if (Tsrf[i, j] + dTs > Tm && Sice[1, i, j] <= eps(Tf))
            Qs = qsat(Ps[i, j], Tm)
            Esrf[i, j] = rho * KWg[i, j] * (Qs - Qa[i, j])
            G[i, j] = Tf(2) * ks1[i, j] * (Tm - Ts1[i, j]) / Ds1[i, j]
            H[i, j] = cp * rho * KH[i, j] * (Tm - Ta[i, j])
            LE[i, j] = Ls * Esrf[i, j]
            Rnet[i, j] = SWsrf[i, j] + trcn[i, j] * LWeff[i, j] - sb * Tm^Tf(4) + (Tf(1) - trcn[i, j]) * sb * Tveg[i, j]^Tf(4)
            dE = Tf(0.0)
            dG = Tf(0.0)
            dH = Tf(0.0)
            dR = Tf(0.0)
            dTs = Tm - Tsrf[i, j]
        end
    end

    # Update surface temperature and fluxes
    Esrf[i, j] = Esrf[i, j] + dE
    G[i, j] = G[i, j] + dG
    H[i, j] = H[i, j] + dH
    LE[i, j] = Lh * Esrf[i, j]
    Rnet[i, j] = Rnet[i, j] + dR
    Tsrf[i, j] = Tsrf[i, j] + dTs

    # Sublimation limited by amount of snow after melt
    Ssub = column_sum(Sice, i, j)
    Ssub -= Melt[i, j] * dt
    if (Ssub > eps(Tf) && Esrf[i, j] * dt > Ssub)
        Esrf[i, j] = Ssub / dt
        LE[i, j] = Ls * Esrf[i, j]
        H[i, j] = Rnet[i, j] - G[i, j] - LE[i, j] - Lf * Melt[i, j]
    end
    Hsrf[i, j] = H[i, j]
    LEsrf[i, j] = LE[i, j]
    Rsrf[i, j] = Rnet[i, j]

    # Define LWsci/LWveg even in open runs
    LWsci[i, j] = LWeff[i, j]
    LWveg[i, j] = Tf(0.0)

    return nothing
end

Base.@propagate_inbounds function energy_balance!(::ForestCover{Tf}, substrate, i, j, state, diag, surface, params, meteo) where {Tf}

    @unpack_constants(Tf)

    (; dt) = params
    (; canh, trcn) = surface
    (; Qcan, Sice, Tcan, Tsrf, Tveg) = state
    (; Esrf, Eveg, G, H, Hsrf, LE, LEsrf, LWsci, LWveg, Melt, Rnet, Rsrf, Ds1, KHa, KHg, KHv, KWg, KWv, ks1, SWsrf, SWveg, Ts1, Tveg0, Qa, LWeff) = diag
    (; Ps, Ta, LW) = meteo

    # Kernel-local scratch variables
    A = zero(MMatrix{4, 4, Tf})
    b = zero(MVector{4, Tf})
    x = zero(MVector{4, Tf})

    # Saturation humidity and density of air
    Qsrf = qsat(Ps[i, j], Tsrf[i, j])
    Lsrf = Ls
    if (Tsrf[i, j] > Tm)
        Lsrf = Lv
    end
    Dsrf = Lsrf * Qsrf / (Rwat * Tsrf[i, j]^Tf(2))
    Qveg = qsat(Ps[i, j], Tveg[i, j])
    Lveg = Ls
    if (Tveg[i, j] > Tm)
        Lveg = Lv
    end
    Dveg = Lveg * Qveg / (Rwat * Tveg[i, j]^Tf(2))
    rho = Ps[i, j] / (Rair * Ta[i, j])

    # Explicit fluxes
    E = rho * KHa[i, j] * (Qcan[i, j] - Qa[i, j])
    Esrf[i, j] = rho * KWg[i, j] * (Qsrf - Qcan[i, j])
    Eveg[i, j] = rho * KWv[i, j] * (Qveg - Qcan[i, j])
    G[i, j] = Tf(2) * ks1[i, j] * (Tsrf[i, j] - Ts1[i, j]) / Ds1[i, j]
    H[i, j] = rho * cp * KHa[i, j] * (Tcan[i, j] - Ta[i, j])
    Hsrf[i, j] = rho * cp * KHg[i, j] * (Tsrf[i, j] - Tcan[i, j])
    Hveg = rho * cp * KHv[i, j] * (Tveg[i, j] - Tcan[i, j])
    LE[i, j] = Lsrf * Esrf[i, j] + Lveg * Eveg[i, j]
    Melt[i, j] = Tf(0)
    Rsrf[i, j] = SWsrf[i, j] + trcn[i, j] * LWeff[i, j] - sb * Tsrf[i, j]^Tf(4) + (Tf(1) - trcn[i, j]) * sb * Tveg[i, j]^Tf(4)        # with near and distant canopy contributions
    Rveg = SWveg[i, j] + (Tf(1) - trcn[i, j]) * (LW[i, j] + sb * Tsrf[i, j]^Tf(4) - Tf(2) * sb * Tveg[i, j]^Tf(4))

    # Surface energy balance increments without melt
    A[1, 1] = Tf(0)
    A[1, 2] = -(KHa[i, j] + KHv[i, j] + KHg[i, j])
    A[1, 3] = KHg[i, j]
    A[1, 4] = KHv[i, j]
    b[1] = (H[i, j] - Hveg - Hsrf[i, j]) / (rho * cp)
    A[2, 1] = -(KHa[i, j] + KWv[i, j] + KWg[i, j])
    A[2, 2] = Tf(0)
    A[2, 3] = Dsrf * KWg[i, j]
    A[2, 4] = Dveg * KWv[i, j]
    b[2] = (E - Eveg[i, j] - Esrf[i, j]) / rho
    A[3, 1] = -Lsrf * rho * KWg[i, j]
    A[3, 2] = -rho * cp * KHg[i, j]
    A[3, 3] = rho * (cp * KHg[i, j] + Lsrf * Dsrf * KWg[i, j]) + Tf(4) * sb * Tsrf[i, j]^Tf(3) + Tf(2) * ks1[i, j] / Ds1[i, j]
    A[3, 4] = -Tf(4) * (Tf(1) - trcn[i, j]) * sb * Tveg[i, j]^Tf(3)
    b[3] = Rsrf[i, j] - Hsrf[i, j] - Lsrf * Esrf[i, j] - G[i, j]
    A[4, 1] = -Lveg * rho * KWv[i, j]
    A[4, 2] = -rho * cp * KHv[i, j]
    A[4, 3] = -Tf(4) * (Tf(1) - trcn[i, j]) * sb * Tsrf[i, j]^Tf(3)
    A[4, 4] = canh[i, j] / dt + rho * (cp * KHv[i, j] + Lveg * Dveg * KWv[i, j]) + Tf(8) * (Tf(1) - trcn[i, j]) * sb * Tveg[i, j]^Tf(3)
    b[4] = Rveg - Hveg - Lveg * Eveg[i, j] - canh[i, j] * (Tveg[i, j] - Tveg0[i, j]) / dt
    ludcmp!(4, A, b, x)
    dQc = x[1]
    dTc = x[2]
    dTs = x[3]
    dTv = x[4]
    dEs = rho * KWg[i, j] * (Dsrf * dTs - dQc)
    dEv = rho * KWv[i, j] * (Dveg * dTv - dQc)
    dGs = Tf(2) * ks1[i, j] * dTs / Ds1[i, j]
    dHs = rho * cp * KHg[i, j] * (dTs - dTc)
    dHv = rho * cp * KHv[i, j] * (dTv - dTc)

    # Surface melting
    if (Tsrf[i, j] + dTs > Tm && Sice[1, i, j] > eps(Sice[1, i, j]))
        Melt[i, j] = column_sum(Sice, i, j) / dt
        b[3] = Rsrf[i, j] - Hsrf[i, j] - Lsrf * Esrf[i, j] - G[i, j] - Lf * Melt[i, j]
        ludcmp!(4, A, b, x)
        dQc = x[1]
        dTc = x[2]
        dTs = x[3]
        dTv = x[4]
        dEs = rho * KWg[i, j] * (Dsrf * dTs - dQc)
        dEv = rho * KWv[i, j] * (Dveg * dTv - dQc)
        dGs = Tf(2) * ks1[i, j] * dTs / Ds1[i, j]
        dHs = rho * cp * KHg[i, j] * (dTs - dTc)
        dHv = rho * cp * KHv[i, j] * (dTv - dTc)
        if (Tsrf[i, j] + dTs < Tm)
            Qsrf = qsat(Ps[i, j], Tm)
            Esrf[i, j] = rho * KWg[i, j] * (Qsrf - Qcan[i, j])
            G[i, j] = Tf(2) * ks1[i, j] * (Tm - Ts1[i, j]) / Ds1[i, j]
            Hsrf[i, j] = rho * cp * KHg[i, j] * (Tm - Tcan[i, j])
            Rsrf[i, j] = SWsrf[i, j] + trcn[i, j] * LWeff[i, j] - sb * Tm^Tf(4) + (Tf(1) - trcn[i, j]) * sb * Tveg[i, j]^Tf(4)
            Rveg = SWveg[i, j] + (Tf(1) - trcn[i, j]) * (LW[i, j] + sb * Tm^Tf(4) - Tf(2) * sb * Tveg[i, j]^Tf(4))
            A[1, 3] = Tf(0)
            b[1] = (H[i, j] - Hveg - Hsrf[i, j]) / (rho * cp)
            A[2, 3] = Tf(0)
            b[2] = (E - Eveg[i, j] - Esrf[i, j]) / rho
            A[3, 3] = Tf(1)
            b[3] = Rsrf[i, j] - Hsrf[i, j] - Lsrf * Esrf[i, j] - G[i, j]
            A[4, 3] = Tf(0)
            b[4] = Rveg - Hveg - Lveg * Eveg[i, j] - canh[i, j] * (Tveg[i, j] - Tveg0[i, j]) / dt
            ludcmp!(4, A, b, x)
            dQc = x[1]
            dTc = x[2]
            Melt[i, j] = x[3] / Lf
            dTv = x[4]
            dTs = Tm - Tsrf[i, j]
            dEs = Tf(0)
            dEv = rho * KWv[i, j] * (Dveg * dTv - dQc)
            dGs = Tf(2) * ks1[i, j] * dTs / Ds1[i, j]
            dHs = Tf(0)
            dHv = rho * cp * KHv[i, j] * (dTv - dTc)
        end
    end

    # Update temperatures and fluxes
    Qcan[i, j] = Qcan[i, j] + dQc
    Tcan[i, j] = Tcan[i, j] + dTc
    Tsrf[i, j] = Tsrf[i, j] + dTs
    Tveg[i, j] = Tveg[i, j] + dTv
    Esrf[i, j] = Esrf[i, j] + dEs
    Eveg[i, j] = Eveg[i, j] + dEv
    Hsrf[i, j] = Hsrf[i, j] + dHs
    Hveg = Hveg + dHv
    G[i, j] = G[i, j] + dGs
    H[i, j] = Hsrf[i, j] + Hveg
    LE[i, j] = Lsrf * Esrf[i, j] + Lveg * Eveg[i, j]
    LEsrf[i, j] = Lsrf * Esrf[i, j]
    Rnet[i, j] = SWsrf[i, j] + SWveg[i, j] + LW[i, j] - trcn[i, j] * sb * Tsrf[i, j]^Tf(4) - (Tf(1) - trcn[i, j]) * sb * Tveg[i, j]^Tf(4)
    LWsci[i, j] = trcn[i, j] * LWeff[i, j] + (Tf(1) - trcn[i, j]) * sb * Tveg[i, j]^Tf(4)
    LWveg[i, j] = (Tf(1) - trcn[i, j]) * (LW[i, j] + sb * Tsrf[i, j]^Tf(4) - Tf(2) * sb * Tveg[i, j]^Tf(4))

    # Sublimation limited by amount of snow after melt
    Ssub = column_sum(Sice, i, j) - Melt[i, j] * dt
    if (Ssub > eps(Tf) && Esrf[i, j] * dt > Ssub)
        Esrf[i, j] = Ssub / dt
        LEsrf[i, j] = Ls * Esrf[i, j]
        Hsrf[i, j] = Rnet[i, j] - G[i, j] - LEsrf[i, j] - Lf * Melt[i, j]
    end

    return nothing
end

# Per-cell delegation for a MixedLandCover: pick the concrete land-cover method from the
# `forestcells` mask. Both branches call concrete, inlined methods (no dynamic dispatch),
# so this stays GPU-safe.

@inline function canopy_snow!(m::MixedLandCover, i, j, state, diag, surface, params)
    if m.forestcells[i, j]
        canopy_snow!(m.forest, i, j, state, diag, surface, params)
    else
        canopy_snow!(m.open, i, j, state, diag, surface, params)
    end
    return nothing
end

@inline function exchange_coefficients!(m::MixedLandCover, i, j, state, diag, surface, params, meteo, z0g)
    if m.forestcells[i, j]
        exchange_coefficients!(m.forest, i, j, state, diag, surface, params, meteo, z0g)
    else
        exchange_coefficients!(m.open, i, j, state, diag, surface, params, meteo, z0g)
    end
    return nothing
end

@inline function solar_radiation!(m::MixedLandCover, i, j, state, diag, surface, meteo)
    if m.forestcells[i, j]
        solar_radiation!(m.forest, i, j, state, diag, surface, meteo)
    else
        solar_radiation!(m.open, i, j, state, diag, surface, meteo)
    end
    return nothing
end

@inline function thermal_radiation!(m::MixedLandCover, i, j, diag, surface, meteo)
    if m.forestcells[i, j]
        thermal_radiation!(m.forest, i, j, diag, surface, meteo)
    else
        thermal_radiation!(m.open, i, j, diag, surface, meteo)
    end
    return nothing
end

Base.@propagate_inbounds function energy_balance!(m::MixedLandCover, substrate::MixedSubstrate, i, j, state, diag, surface, params, meteo)
    if m.forestcells[i, j]
        energy_balance!(m.forest, substrate.soil, i, j, state, diag, surface, params, meteo)
    elseif substrate.icecells[i, j]
        energy_balance!(m.open, substrate.ice, i, j, state, diag, surface, params, meteo)
    else
        energy_balance!(m.open, substrate.soil, i, j, state, diag, surface, params, meteo)
    end
    return nothing
end

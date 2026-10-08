@kwdef struct DiagnosticAlbedo{Tf} <: AbstractAlbedo{Tf}
    amin::Tf = 0.6                            # Minimum albedo for melting snow (-)
    amax::Tf = 0.86                           # Maximum albedo for fresh snow (-)
    Talb::Tf = -2                             # Albedo decay temperature threshold (C)
end

@kwdef struct DecayAlbedo{Tf, GT, MF <: AbstractMatrix{<:AbstractFloat}} <: AbstractAlbedo{Tf}
    grid::GT
    amin::Tf = 0.6                            # Minimum albedo for melting snow (-)
    tcld::Tf = 3600 * 1000                    # Cold snow albedo decay time scale (s)
    tmlt::Tf = 3600 * 100                     # Melting snow albedo decay time scale (s)
    adfs::Tf = 3                              # Albedo adjustment, shortwave (-)
    adfl::Tf = 2                              # Albedo adjustment, longwave (-)
    Sfmin::Tf = 10                            # Minimum snowfall over 24h to refresh albedo (kg/m^2)
    afs::MF = 0.86 * ones(grid.Nx, grid.Ny)   # Maximum albedo for fresh snow (-)
end

@kwdef struct PrognosticAlbedo{Tf, GT, MF <: AbstractMatrix{<:AbstractFloat}} <: AbstractAlbedo{Tf}
    grid::GT
    ALRADT::Bool = true                       # Aspect-dependent decay tuning
    adm::Tf = 100                             # Melting snow albedo decay time (h)
    amin::Tf = 0.6                            # Minimum albedo for melting snow (-)
    Sfmin::Tf = 10                            # Minimum snowfall over 24h to refresh albedo (kg/m^2)
    afs::MF = 0.86 * ones(grid.Nx, grid.Ny)   # Maximum albedo for fresh snow (-)
    adc::MF = 1000 * ones(grid.Nx, grid.Ny)   # Cold snow albedo decay time (h)
end

DiagnosticAlbedo{Tf}(grid::Grid; kwargs...) where {Tf} = DiagnosticAlbedo{Tf}(; kwargs...)
DecayAlbedo{Tf}(grid::Grid; kwargs...) where {Tf} = DecayAlbedo{Tf, typeof(grid), Matrix{Tf}}(; grid, kwargs...)
PrognosticAlbedo{Tf}(grid::Grid; kwargs...) where {Tf} = PrognosticAlbedo{Tf, typeof(grid), Matrix{Tf}}(; grid, kwargs...)

@adapt_structure DecayAlbedo
@adapt_structure PrognosticAlbedo

"""
    snow_albedo!(scheme, i, j, state, surface, meteo, params, summer_decay)

Update snow albedo `albs[i, j]` for cell `(i, j)` implemented for every 
`AbstractAlbedo`.
"""
function snow_albedo! end

@inline function snow_albedo!(scheme::DiagnosticAlbedo{Tf}, i, j, state, surface, meteo, params, summer_decay) where {Tf}
    @unpack_constants(Tf)
    (; albs, Tsrf) = state
    afs_loc = scheme.amax
    a = scheme.amin + (afs_loc - scheme.amin) * (Tsrf[i, j] - Tm) / scheme.Talb
    a = max(a, min(afs_loc, scheme.amin))
    a = min(a, max(afs_loc, scheme.amin))
    albs[i, j] = a
    return nothing
end

@inline function snow_albedo!(scheme::DecayAlbedo{Tf}, i, j, state, surface, meteo, params, summer_decay) where {Tf}
    @unpack_constants(Tf)
    (; albs, Tsrf) = state
    (; fveg, trcn, fsky) = surface
    (; Sdir, Sdif, Sf, Tv) = meteo
    (; dt) = params
    afs_loc = scheme.afs[i, j]

    tau = scheme.tcld
    if (Tsrf[i, j] >= Tm)
        tau = scheme.tmlt
    end

    # Melt-season decay is a fixed 70 h, overriding both scheme timescales
    if summer_decay
        tau = Tf(70.0) * Tf(3600.0)
    end

    if fveg[i, j] > Tf(0) && Sdir[i, j] > eps(Tf)
        tau = tau / ((Tf(1) - trcn[i, j] * fsky[i, j]) * (Tf(1) + scheme.adfl * Tv[i, j]) + scheme.adfs * Tv[i, j])
    elseif fveg[i, j] > Tf(0) && Sdif[i, j] > eps(Tf)
        tau = tau / ((Tf(1) - trcn[i, j] * fsky[i, j]) + scheme.adfs * trcn[i, j] * fsky[i, j])
    elseif (fveg[i, j] > Tf(0) && (Sdir[i, j] + Sdif[i, j] <= eps(Tf)))
        tau = tau / (Tf(2.0) - trcn[i, j] * fsky[i, j])
    end

    rt = Tf(1) / tau + Sf[i, j] / scheme.Sfmin
    alim = (scheme.amin / tau + Sf[i, j] * afs_loc / scheme.Sfmin) / rt
    a = alim + (albs[i, j] - alim) * exp(-rt * dt)
    if (a < min(afs_loc, scheme.amin))
        a = min(afs_loc, scheme.amin)
    end
    if (a > max(afs_loc, scheme.amin))
        a = max(afs_loc, scheme.amin)
    end
    albs[i, j] = a
    return nothing
end

@inline function snow_albedo!(scheme::PrognosticAlbedo{Tf}, i, j, state, surface, meteo, params, summer_decay) where {Tf}
    @unpack_constants(Tf)
    (; albs, Tsrf, Sice, Sliq) = state
    (; Sdir, Sdird, Sf, Sf24h) = meteo
    (; dt) = params
    adc_loc = scheme.adc[i, j]
    adm_loc = scheme.adm
    afs_loc = scheme.afs[i, j]

    SWEtmp = zero(Tf)
    for si in 1:size(Sice, 1)
        SWEtmp += Sice[si, i, j] + Sliq[si, i, j]
    end

    # Aspect-dependent albedo tuning
    if scheme.ALRADT
        if ((Sdir[i, j] > eps(Tf)) && (Sdird[i, j] < Sdir[i, j]))
            adm_loc = adm_loc * (Sdird[i, j]) / (Sdir[i, j])
            adc_loc = adc_loc * (Sdird[i, j]) / (Sdir[i, j])
            if (adm_loc < eps(Tf))
                adm_loc = eps(Tf)
            end
            if (adc_loc < eps(Tf))
                adc_loc = eps(Tf)
            end
        end
    end

    # Temperature dependent albedo update
    a = albs[i, j]
    if (Tsrf[i, j] >= Tm)
        a = (a - scheme.amin) * exp(-(dt / Tf(3600)) / adm_loc) + scheme.amin
    else
        a = a - (dt / Tf(3600)) / adc_loc
    end

    # Reduce albedo for thin and patchy snow cover
    if (SWEtmp < Tf(75.0))
        afs_loc *= Tf(0.8)
    end

    # Reset to fresh snow albedo
    if ((Sf[i, j] * dt) > Tf(0.0) && Sf24h[i, j] > scheme.Sfmin)
        a = afs_loc
    else
        a = a + (afs_loc - a) * Sf[i, j] * dt / scheme.Sfmin
    end

    if (a > afs_loc)
        a = afs_loc
    end
    if (a < scheme.amin)
        a = scheme.amin
    end
    albs[i, j] = a
    return nothing
end

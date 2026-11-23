@kwdef struct AgeCompaction{Tf} <: AbstractCompaction{Tf}
    rmlt::Tf = 500              # Maximum density for melting snow (kg/m^3)
    rcld::Tf = 300              # Maximum density for cold snow (kg/m^3)
    trho::Tf = 3600 * 200       # Snow compaction time scale (s)
end

@kwdef struct OverburdenCompaction{Tf} <: AbstractCompaction{Tf}
    eta0::Tf = 3.7e7            # Reference snow viscosity (Pa s)
    snda::Tf = 2.8e-6           # Thermal metamorphism parameter (1/s)
    rhos_max::Tf = 750          # Maximum snow density (kg/m^3)
end

"""
$(TYPEDEF)

Crocus (B92) overburden compaction with an optional wind compaction of surface snow
(Groot Zwaaftink et al., 2013).

# Fields

$(TYPEDFIELDS)
"""
@kwdef struct CrocusCompaction{Tf} <: AbstractCompaction{Tf}
    "Reference snow viscosity (Pa s)"
    eta1::Tf = 7.62237e6
    "Temperature factor (1/K)"
    a_eta::Tf = 0.1
    "First density factor (m^3/kg)"
    b_eta::Tf = 0.023
    "Second density factor (kg/m^3)"
    c_eta::Tf = 250
    "Maximum snow density (kg/m^3)"
    rhos_max::Tf = 750
    "Liquid-water viscosity factor (-)"
    f_liq::Tf = 600
    "Wind compaction mode: 0 off, 1 cubic, 2 linear in the excess wind speed"
    wndens::Int = 0
    "Threshold wind speed at 3 m for wind compaction (m/s)"
    U_thresh::Tf = 5.0
    "Wind compaction enhancement parameter (s^3/m^3)"
    A_enh0::Tf = 5.0
    "Maximum depth affected by wind compaction (m)"
    d_max::Tf = 0.07
end

AgeCompaction{Tf}(grid::Grid; kwargs...) where {Tf} = AgeCompaction{Tf}(; kwargs...)
OverburdenCompaction{Tf}(grid::Grid; kwargs...) where {Tf} = OverburdenCompaction{Tf}(; kwargs...)
CrocusCompaction{Tf}(grid::Grid; kwargs...) where {Tf} = CrocusCompaction{Tf}(; kwargs...)

"""
    compact_snow!(scheme, i, j, state, diag, surface, params)

Compact the snow column at cell `(i, j)`: rescale the layer thicknesses `Ds` in
place to the compacted density, for every layer. `diag`/`surface` carry the
effective wind `Uaeff` and roughness `z0_snow` used by the `CrocusCompaction`
wind-compaction option; the other schemes ignore them.
"""
function compact_snow! end

@inline function compact_snow!(scheme::AgeCompaction{Tf}, i, j, state, diag, surface, params) where {Tf}
    @unpack_constants(Tf)
    (; Ds, Sice, Sliq, fsnow, Nsnow, Tsnow) = state
    (; dt) = params
    (; rmlt, rcld, trho) = scheme
    for k in 1:Nsnow[i, j]
        if (Ds[k, i, j] > eps(Tf))
            rhos = (Sice[k, i, j] + Sliq[k, i, j]) / Ds[k, i, j] / fsnow[i, j]
            if (Tsnow[k, i, j] >= Tm)
                if (rhos < rmlt)
                    rhos = rmlt + (rhos - rmlt) * exp(-dt / trho)
                end
            else
                if (rhos < rcld)
                    rhos = rcld + (rhos - rcld) * exp(-dt / trho)
                end
            end
            Ds[k, i, j] = (Sice[k, i, j] + Sliq[k, i, j]) / rhos / fsnow[i, j]
        end
    end
    return nothing
end

@inline function compact_snow!(scheme::OverburdenCompaction{Tf}, i, j, state, diag, surface, params) where {Tf}
    @unpack_constants(Tf)
    (; Ds, Sice, Sliq, fsnow, Nsnow, Tsnow) = state
    (; dt) = params
    (; eta0, snda, rhos_max) = scheme
    mass = Tf(0.0)
    for k in 1:Nsnow[i, j]
        mass = mass + Tf(0.5) * (Sice[k, i, j] + Sliq[k, i, j]) / fsnow[i, j]
        if (Ds[k, i, j] > eps(Tf))
            rhos = (Sice[k, i, j] + Sliq[k, i, j]) / Ds[k, i, j] / fsnow[i, j]
            rhos = rhos + (rhos * grav * mass * dt / (eta0 * exp(-(Tsnow[k, i, j] - Tm) / Tf(12.4) + rhos / Tf(55.6))) + dt * rhos * snda * exp((Tsnow[k, i, j] - Tm) / Tf(23.8) - max(rhos - Tf(150), Tf(0.0)) / Tf(21.7)))
            rhos = min(rhos, rhos_max)
            Ds[k, i, j] = (Sice[k, i, j] + Sliq[k, i, j]) / rhos / fsnow[i, j]
        end
        mass = mass + Tf(0.5) * (Sice[k, i, j] + Sliq[k, i, j]) / fsnow[i, j]
    end
    return nothing
end

@inline function compact_snow!(scheme::CrocusCompaction{Tf}, i, j, state, diag, surface, params) where {Tf}
    @unpack_constants(Tf)
    (; Ds, Sice, Sliq, fsnow, Nsnow, Tsnow) = state
    (; dt, zU) = params
    (; eta1, a_eta, b_eta, c_eta, rhos_max, f_liq, wndens, U_thresh, A_enh0, d_max) = scheme
    (; Uaeff) = diag
    (; z0_snow) = surface
    mass = Tf(0.0)
    d_bottom = Tf(0.0)
    # Wind speed at 3 m, where the wind compaction parameterization is valid
    U_new = wndens == 0 ? Tf(0) : Uaeff[i, j] * log(Tf(3) / z0_snow[i, j]) / log(zU / z0_snow[i, j])
    for k in 1:Nsnow[i, j]
        mass = mass + Tf(0.5) * (Sice[k, i, j] + Sliq[k, i, j]) / fsnow[i, j]
        d_bottom = d_bottom + Ds[k, i, j]
        d_mid = d_bottom - Ds[k, i, j] / Tf(2)
        if (Ds[k, i, j] > eps(Tf))
            rhos = (Sice[k, i, j] + Sliq[k, i, j]) / Ds[k, i, j] / fsnow[i, j]
            f1 = Tf(1) / (Tf(1) + f_liq * Sliq[k, i, j] / (rho_wat * Ds[k, i, j] * fsnow[i, j]))
            f2 = Tf(1.0)
            eta = f1 * f2 * eta1 * (rhos / c_eta) * exp(a_eta * (Tm - Tsnow[k, i, j]) + b_eta * rhos)
            if wndens == 0
                rhos = rhos + rhos * grav * mass * dt / eta
            else
                # Wind compaction of surface snow, Groot Zwaaftink et al. (2013), eq. 2-3
                g_wind = Tf(0)
                if U_new > U_thresh && d_mid <= d_max
                    if wndens == 1
                        A_enh = Tf(2.7) * A_enh0 * (Tf(1) - d_mid / (Tf(1.25) * d_max))
                        g_wind = A_enh * (U_new - U_thresh)^3
                    elseif wndens == 2
                        g_wind = A_enh0 * (U_new - U_thresh)
                    end
                end
                rhos = rhos + rhos * dt * (Tf(1) + g_wind) * grav * mass / eta
            end
            rhos = min(rhos, rhos_max)
            Ds[k, i, j] = (Sice[k, i, j] + Sliq[k, i, j]) / rhos / fsnow[i, j]
        end
        mass = mass + Tf(0.5) * (Sice[k, i, j] + Sliq[k, i, j]) / fsnow[i, j]
    end
    return nothing
end

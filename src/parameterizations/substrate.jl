@kwdef struct SoilSubstrate{Tf} <: AbstractSubstrate{Tf}
    gsat::Tf = 0.01             # Surface conductance for saturated soil (m/s)
end

@kwdef struct IceSubstrate{Tf} <: AbstractSubstrate{Tf}
    gsat::Tf = 0.01             # Surface conductance for glacier ice (m/s)
end

SoilSubstrate{Tf}(grid::Grid; kwargs...) where {Tf} = SoilSubstrate{Tf}(; kwargs...)
IceSubstrate{Tf}(grid::Grid; kwargs...) where {Tf} = IceSubstrate{Tf}(; kwargs...)

# Per-cell mix of soil and glacier-ice substrate for a single run. `icecells` is a per-cell mask:
# physics uses the ice scheme where it is true and the soil scheme elsewhere. Carries a grid-sized
# array like the per-cell albedo schemes, so it is `@adapt_structure`d.
struct MixedSubstrate{Tf, S, I, MB} <: AbstractSubstrate{Tf}
    soil::S
    ice::I
    icecells::MB
end
MixedSubstrate(soil::SoilSubstrate{Tf}, ice::IceSubstrate{Tf}, icecells::AbstractMatrix{Bool}) where {Tf} =
    MixedSubstrate{Tf, typeof(soil), typeof(ice), typeof(icecells)}(soil, ice, icecells)
@adapt_structure MixedSubstrate

"""
    soil_properties!(substrate, i, j, state, diag, surface, grid, params)

Fill the soil heat capacity `diag.csoil[1:Nsoil, i, j]`, thermal conductivity
`diag.ksoil[1:Nsoil, i, j]` and surface moisture conductance `diag.gs1[i, j]` for cell
`(i, j)`. Called from the `thermal!` kernel.
"""
function soil_properties! end

@inline function soil_properties!(s::IceSubstrate{Tf}, i, j, state, diag, surface, grid, params) where {Tf}
    @unpack_constants(Tf)
    (; Dzsoil, Nsoil) = grid
    (; gsat) = s
    (; csoil, ksoil, gs1) = diag

    for k in 1:Nsoil
        # hcap_ice is a specific heat capacity and needs converting to a volumetric value
        csoil[k, i, j] = hcap_ice * rho_ice * Dzsoil[k]
        ksoil[k, i, j] = hcon_ice
        # An ice surface behaves like saturated soil for surface moisture conductance
        gs1[i, j] = gsat
    end
    return nothing
end

@inline function soil_properties!(s::SoilSubstrate{Tf}, i, j, state, diag, surface, grid, params) where {Tf}
    @unpack_constants(Tf)
    (; Dzsoil, Nsoil) = grid
    (; gsat) = s
    (; b, hcap_soil, hcon_soil, sathh, Vcrit, Vsat) = surface
    (; theta, Tsoil) = state
    (; csoil, ksoil, gs1) = diag

    dPsidT = -rho_ice * Lf / (rho_wat * grav * Tm)

    for k in 1:Nsoil
        csoil[k, i, j] = hcap_soil[i, j] * Dzsoil[k]
        ksoil[k, i, j] = hcon_soil[i, j]
        if (theta[k, i, j] > eps(Tf))
            dthudT = Tf(0.0)
            sthu = theta[k, i, j]
            sthf = Tf(0.0)
            Tc = Tsoil[k, i, j] - Tm
            Tmax = Tm + (sathh[i, j] / dPsidT) * (Vsat[i, j] / theta[k, i, j])^b[i, j]
            if (Tsoil[k, i, j] < Tmax)
                dthudT = (-dPsidT * Vsat[i, j] / (b[i, j] * sathh[i, j])) * (dPsidT * Tc / sathh[i, j])^(Tf(-1) / b[i, j] - Tf(1))
                sthu = Vsat[i, j] * (dPsidT * Tc / sathh[i, j])^(Tf(-1) / b[i, j])
                sthu = min(sthu, theta[k, i, j])
                sthf = (theta[k, i, j] - sthu) * rho_wat / rho_ice
            end
            Mf = rho_ice * Dzsoil[k] * sthf
            Mu = rho_wat * Dzsoil[k] * sthu
            csoil[k, i, j] = hcap_soil[i, j] * Dzsoil[k] + hcap_ice * Mf + hcap_wat * Mu + rho_wat * Dzsoil[k] * ((hcap_wat - hcap_ice) * Tc + Lf) * dthudT
            Smf = rho_ice * sthf / (rho_wat * Vsat[i, j])
            Smu = sthu / Vsat[i, j]
            thice = Tf(0.0)
            if (Smf > eps(Tf))
                thice = Vsat[i, j] * Smf / (Smu + Smf)
            end
            thwat = Tf(0.0)
            if (Smu > eps(Tf))
                thwat = Vsat[i, j] * Smu / (Smu + Smf)
            end
            hcon_sat = hcon_soil[i, j] * (hcon_wat^thwat) * (hcon_ice^thice) / (hcon_air^Vsat[i, j])
            ksoil[k, i, j] = (hcon_sat - hcon_soil[i, j]) * (Smf + Smu) + hcon_soil[i, j]
            if (k == 1)
                gs1[i, j] = gsat * max((Smu * Vsat[i, j] / Vcrit[i, j])^Tf(2), Tf(1.0))
            end

        end

    end
    return nothing
end

@inline function soil_properties!(m::MixedSubstrate, i, j, state, diag, surface, grid, params)
    if m.icecells[i, j]
        soil_properties!(m.ice, i, j, state, diag, surface, grid, params)
    else
        soil_properties!(m.soil, i, j, state, diag, surface, grid, params)
    end
    return nothing
end

"""
    cap_soil_temperature!(substrate, i, j, state, grid)

Cap the substrate temperature at cell `(i, j)` such that glacier ice cannot exceed
 the melting point, so `IceSubstrate` clamps it there and discards the excess energy.
"""
function cap_soil_temperature! end

@inline cap_soil_temperature!(::SoilSubstrate, i, j, state, grid) = nothing

@inline function cap_soil_temperature!(::IceSubstrate{Tf}, i, j, state, grid) where {Tf}
    @unpack_constants(Tf)
    (; Nsoil) = grid
    (; Tsoil) = state

    for k in 1:Nsoil
        Tsoil[k, i, j] = min(Tsoil[k, i, j], Tm)
    end
    return nothing
end

@inline function cap_soil_temperature!(m::MixedSubstrate, i, j, state, grid)
    if m.icecells[i, j]
        cap_soil_temperature!(m.ice, i, j, state, grid)
    else
        cap_soil_temperature!(m.soil, i, j, state, grid)
    end
    return nothing
end

"""
    cap_initial_temperatures!(substrate, state, surface, grid)

Cap the initial surface and soil temperatures at the melting point on glacier ice: everywhere for an
`IceSubstrate`, on the `icecells` for a `MixedSubstrate`, and nowhere for a `SoilSubstrate`. Called
once from `build_state`.
"""
function cap_initial_temperatures! end

cap_initial_temperatures!(::SoilSubstrate, state, surface, grid) = nothing

function cap_initial_temperatures!(::IceSubstrate{Tf}, state, surface, grid) where {Tf}
    Tm = get_constants(Tf).Tm
    state.Tsrf .= min.(state.Tsrf, Tm)
    state.Tsoil .= min.(state.Tsoil, Tm)
    return nothing
end

function cap_initial_temperatures!(m::MixedSubstrate, state, surface, grid)
    Tm = get_constants(eltype(state.Tsrf)).Tm
    glacier = m.icecells
    state.Tsrf[glacier] .= min.(state.Tsrf[glacier], Tm)
    for k in 1:grid.Nsoil
        Tsoilk = @view state.Tsoil[k, :, :]
        Tsoilk[glacier] .= min.(Tsoilk[glacier], Tm)
    end
    return nothing
end

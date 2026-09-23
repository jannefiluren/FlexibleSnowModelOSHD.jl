"""
$(TYPEDSIGNATURES)

Snow albedo calculations, surface and canopy net shortwave radiation, terrain correction of
longwave radiation for open terrain.

# Arguments
- `fsm::FSM`: Model state structure
- `meteo::MET`: Current meteorological conditions
- `t`: Current simulation time
"""
function radiation!(fsm::FSM{Tf}, meteo::MET{Tf}, t) where {Tf <: Real}

    (; land_cover, snow_albedo) = fsm.physics

    summer_decay = Dates.value(Month(t)) > 4 && Dates.value(Month(t)) < 10

    backend = get_backend(fsm.state.albs)
    kernel! = radiation_kernel!(backend)
    kernel!(
        fsm.state, fsm.diag, fsm.surface, fsm.params, meteo,
        land_cover, snow_albedo, summer_decay;
        ndrange = (Int(fsm.grid.Nx), Int(fsm.grid.Ny))
    )
    KernelAbstractions.synchronize(backend)

    return nothing
end

@kernel function radiation_kernel!(
        state, diag, surface, params::Parameters{Tf}, meteo,
        land_cover::AbstractLandCover{Tf},
        albedo::AbstractAlbedo{Tf}, summer_decay::Bool,
    ) where {Tf}

    i, j = @index(Global, NTuple)

    (; active, alb0) = surface
    (; albs, fsnow) = state

    if active[i, j]

        snow_albedo!(albedo, i, j, state, surface, meteo, params, summer_decay)

        # Bare ground shows through once the snow has gone
        if (fsnow[i, j] <= eps(Tf))
            albs[i, j] = alb0[i, j]
        end

        # Surface albedo, shortwave transmission and thermal emission from surroundings
        solar_radiation!(land_cover, i, j, state, diag, surface, meteo)
        thermal_radiation!(land_cover, i, j, diag, surface, meteo)

    end
end

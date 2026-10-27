"""
$(TYPEDSIGNATURES)

Eddy diffusivities for turbulent transfer of heat and moisture between the ground,
the canopy and the atmosphere.

# Arguments
- `fsm::FSM`: Model state structure
- `meteo::MET`: Current meteorological conditions
"""
function surface_exchange_coefficients!(fsm::FSM{Tf}, meteo::MET{Tf}) where {Tf <: Real}

    (; land_cover, snow_fraction, snow_roughness) = fsm.physics

    backend = get_backend(fsm.state.Tsrf)
    kernel! = surface_exchange_coefficients_kernel!(backend)
    kernel!(
        fsm.state, fsm.diag, fsm.surface, fsm.params, meteo,
        land_cover, snow_fraction, snow_roughness;
        ndrange = (Int(fsm.grid.Nx), Int(fsm.grid.Ny))
    )
    KernelAbstractions.synchronize(backend)

    return nothing
end

@kernel function surface_exchange_coefficients_kernel!(
        state, diag, surface, params::Parameters{Tf}, meteo,
        land_cover::AbstractLandCover{Tf},
        snow_fraction::AbstractSnowFraction{Tf},
        snow_roughness::AbstractSnowRoughness{Tf},
    ) where {Tf}

    i, j = @index(Global, NTuple)

    (; active) = surface

    if active[i, j]

        z0g = ground_roughness(snow_fraction, snow_roughness, i, j, state, surface)

        exchange_coefficients!(land_cover, i, j, state, diag, surface, params, meteo, z0g)

    end
end

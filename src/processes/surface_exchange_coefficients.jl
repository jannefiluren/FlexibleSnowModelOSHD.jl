# Computation of turbulent eddy diffusivities for heat and moisture

"""
    surface_exchange_coefficients!(fsm, meteo)

Eddy diffusivities for turbulent transfer of heat and moisture between the ground,
the canopy and the atmosphere.

# Arguments
- `fsm::FSM`: Model state structure
- `meteo::MET`: Current meteorological conditions
"""
function surface_exchange_coefficients!(fsm::FSM{Tf}, meteo::MET{Tf}) where {Tf <: Real}

    (; surface_layer, snow_fraction) = fsm.physics

    backend = get_backend(fsm.state.Tsrf)
    kernel! = surface_exchange_coefficients_kernel!(backend)
    kernel!(
        fsm.state, fsm.diag, fsm.surface, fsm.params, meteo,
        surface_layer, snow_fraction;
        ndrange = (Int(fsm.grid.Nx), Int(fsm.grid.Ny))
    )
    KernelAbstractions.synchronize(backend)

    return nothing
end

@kernel function surface_exchange_coefficients_kernel!(
        state, diag, surface, params::Parameters{Tf}, meteo,
        surface_layer::AbstractSurfaceLayer{Tf},
        snow_fraction::AbstractSnowFraction{Tf},
    ) where {Tf}

    i, j = @index(Global, NTuple)

    (; tthresh) = params
    (; tilefrac) = surface

    if (tilefrac[i, j] >= tthresh)

        z0g = ground_roughness(snow_fraction, i, j, state, surface)

        exchange_coefficients!(surface_layer, i, j, state, diag, surface, params, meteo, z0g)

    end
end

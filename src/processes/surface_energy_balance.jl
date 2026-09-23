# Surface energy balance

"""
$(TYPEDSIGNATURES)

Surface energy balance solution, coupled to the canopy where there is one.

# Arguments
- `fsm::FSM`: Model state structure
- `meteo::MET`: Current meteorological conditions
"""
function surface_energy_balance!(fsm::FSM{Tf}, meteo::MET{Tf}) where {Tf <: Real}

    (; land_cover, substrate) = fsm.physics

    backend = get_backend(fsm.state.Tsrf)
    kernel! = surface_energy_balance_kernel!(backend)
    kernel!(
        fsm.state, fsm.diag, fsm.surface, fsm.params, meteo,
        land_cover, substrate;
        ndrange = (Int(fsm.grid.Nx), Int(fsm.grid.Ny))
    )
    KernelAbstractions.synchronize(backend)

    return nothing
end

# inbounds = true (not a raw @inbounds block, which miscompiles the KA CPU kernel) keeps the forest solver's scratch off the heap
@kernel inbounds = true function surface_energy_balance_kernel!(
        state, diag, surface, params::Parameters{Tf}, meteo,
        land_cover::AbstractLandCover{Tf}, substrate::AbstractSubstrate{Tf},
    ) where {Tf}

    i, j = @index(Global, NTuple)

    (; active) = surface

    if active[i, j]

        energy_balance!(land_cover, substrate, i, j, state, diag, surface, params, meteo)

    end
end

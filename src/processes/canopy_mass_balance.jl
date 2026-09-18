# Canopy process: snow interception, sublimation and unloading

"""
    canopy_mass_balance!(fsm)

Snow interception, sublimation, and unloading from the vegetation canopy.

# Arguments
- `fsm::FSM`: Model state structure
"""
function canopy_mass_balance!(fsm::FSM{Tf}) where {Tf <: Real}

    (; land_cover) = fsm.physics

    backend = get_backend(fsm.state.Sveg)
    kernel! = canopy_mass_balance_kernel!(backend)
    kernel!(
        fsm.state, fsm.diag, fsm.surface, fsm.params,
        land_cover;
        ndrange = (Int(fsm.grid.Nx), Int(fsm.grid.Ny))
    )
    KernelAbstractions.synchronize(backend)

    return nothing
end

@kernel function canopy_mass_balance_kernel!(
        state, diag, surface, params::Parameters{Tf},
        land_cover::AbstractLandCover{Tf},
    ) where {Tf}

    i, j = @index(Global, NTuple)

    (; active) = surface

    if active[i, j]

        canopy_snow!(land_cover, i, j, state, diag, surface, params)

    end
end

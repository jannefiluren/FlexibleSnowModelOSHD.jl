"""
$(TYPEDEF)

Top-level model container bundling the grid, parameters, and the surface, state and
diagnostic fields together with the selected physics schemes.

# Fields

$(TYPEDFIELDS)
"""
mutable struct FSM{Tf, G, P, L, S, D, PH}
    "Grid definition"
    grid::G
    "Scalar model parameters"
    params::P
    "Static surface, terrain, canopy and soil properties"
    surface::L
    "Prognostic model state"
    state::S
    "Per-cell diagnostic fields"
    diag::D
    "Selected physics schemes"
    physics::PH
end

# Positional constructor from already-built sub-structs
function FSM(
        grid::Grid, params::Parameters{Tf}, surface::Surface, state::State,
        diag::Diagnostics, physics
    ) where {Tf}
    return FSM{
        Tf, typeof(grid), typeof(params), typeof(surface), typeof(state),
        typeof(diag), typeof(physics),
    }(grid, params, surface, state, diag, physics)
end

function Base.show(io::IO, fsm::FSM{Tf}) where {Tf}
    print(io, "FSM", '\n')
    print(io, "├── Precision: ", Tf, '\n')
    print(io, "├── Nx: ", fsm.grid.Nx, ", Ny: ", fsm.grid.Ny, '\n')
    print(io, "└── Active cells: ", sum(fsm.surface.active), '\n')
    print(io, "Physics options:", '\n')
    for (i, (name, p)) in enumerate(pairs(fsm.physics))
        start = i == length(fsm.physics) ? "└── " : "├── "
        print(io, "   $start", name, ": ", nameof(typeof(p)), '\n')
    end
    return nothing
end

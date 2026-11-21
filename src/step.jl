"""
    step!(fsm, met, t; transport = nothing)

Execute one complete physics time step of the snow model.

This function encapsulates the standard model execution sequence:
1. Meteorological data processing
2. Radiation calculations
3. Thermal property updates
4. Iterative surface energy balance
5. Canopy mass-balance processes
6. Snow column processes
7. Soil thermal processes
8. Horizontal snow transport

# Arguments
- `fsm::FSM`: Model state structure
- `met::MET`: Current meteorological conditions
- `t::DateTime`: Current simulation time
- `transport::Union{SnowTransport, Nothing}` (keyword): when a workspace is passed, run the
  horizontal snow-transport step ([`transport!`](@ref)) at the end of the step, which then
  also rolls the 14-day snow-cover history; when `nothing` (default), the step is identical
  to a run without transport.
"""
function step!(fsm::FSM{Tf}, met::MET{Tf}, t; transport = nothing) where {Tf}

    # 1. Meteorological data processing
    drive!(fsm, met)

    # 2. Radiation calculations
    radiation!(fsm, met, t)

    # 3. Thermal property updates
    thermal!(fsm)

    # 4. Iterative energy balance solution
    for _ in 1:fsm.params.Nitr
        surface_exchange_coefficients!(fsm, met)
        surface_energy_balance!(fsm, met)
    end

    # 5. Canopy interception / unloading
    canopy_mass_balance!(fsm)

    # 6. Snow processes (transport, if active, rolls the snow-cover history instead)
    snow!(fsm, met, t; update_hist = !transport_active(transport) && history_roll_hour(t))

    # 7. Soil thermal processes
    soil!(fsm)

    # 8. Horizontal snow transport
    transport === nothing || transport!(fsm, met, transport, t)

    return nothing
end

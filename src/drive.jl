"""
    drive!(fsm, meteo)

Meteorological data preprocessing and unit conversions.

# Arguments
- `fsm::FSM`: Model state structure
- `meteo::MET`: Current meteorological conditions (modified in-place)
"""
function drive!(fsm::FSM, meteo::MET{Tf,Ti}) where {Tf<:Real,Ti<:Integer}

  @unpack_constants(Tf)

  @unpack dt = fsm

  @unpack Tc, es, Qa, Ua, Sf, Rf, Ta, RH, Ps = meteo

  Ua .= Tf(0.7) .* Ua
  Ua .= max.(Ua, Tf(0.1))

  Sf .= Sf ./ dt
  Rf .= Rf ./ dt
  Tc .= Ta .- Tm
  es .= e0 .* exp.(Tf(17.5043) .* Tc ./ (Tf(241.3) .+ Tc))
  Qa .= (RH ./ 100) .* eps_fsm .* es ./ Ps

end

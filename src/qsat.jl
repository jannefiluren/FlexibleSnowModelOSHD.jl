function qsat(P::Tf, T::Tf) where {Tf <: Real}

  eps_fsm = Tf(0.622)   # Ratio of molecular weights of water and dry air
  Tm = Tf(273.15)       # Melting point (K)
  e0 = Tf(610.78)       # Saturation vapour pressure at Tm (Pa)

  Tc = T - Tm
  if (Tc > Tf(0))
    es = e0 * exp(Tf(17.5043) * Tc / (Tf(241.3) + Tc))
  else
    es = e0 * exp(Tf(22.4422) * Tc / (Tf(272.186) + Tc))
  end
  Qs = eps_fsm * es / P

  return Qs

end
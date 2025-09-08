function fresh_snow_density!(fsm::FSM{Tf, Ti}, Ta, Ua, dem) where {Tf<:Real, Ti<:Integer}
  
  @unpack_constants(Tf)

  @unpack rho0, rhob, rhoc, rhof, rhos_min = fsm

  @unpack DENSTY, OSHDTN = fsm

  if (DENSTY == 0)
    rhonew = rho0
  else
    if (OSHDTN == 0)
      # Initial formulation
      rhonew = max(rhof + rhob * (Ta - Tm) + rhoc * Ua^Tf(0.5), Tf(50.0))
    else # OSHDTN == 1
      # New formulation with decompaction
      rhonew = rhof + rhob * (Ta - Tm) + rhoc * Ua^Tf(0.5)
      if (dem <= Tf(1000))
        t_decompaction = Tf(24.0)
      elseif (dem > Tf(4000))
        t_decompaction = Tf(0.0)
      else
        t_decompaction = Tf(24) + (dem - Tf(1000)) / (Tf(4000) - Tf(1000)) * (Tf(0) - Tf(24))
      end
      rhonew = Tf(300) + (rhonew - Tf(300)) * exp(t_decompaction / Tf(100))
      rhonew = max(rhonew, rhos_min)
    end
  end
  
  return rhonew
  
end
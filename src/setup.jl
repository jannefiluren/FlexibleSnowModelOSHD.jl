function setup(Tf, Ti, landuse::Dict, Nx::Int, Ny::Int; SNFRAC=nothing, TILE="open")

  # Create fsm object
  fsm = FSM{Tf,Ti}(Nx=Nx, Ny=Ny)

  # Model configuration overrides
  if SNFRAC !== nothing
    fsm.SNFRAC = SNFRAC
  end

  # Forest tiles settings
  if TILE == "forest"
    fsm.CANMOD = 1
    fsm.EXCHNG = 2
    fsm.SNFRAC = 4
    fsm.ZOFFST = 1
  end

  # Modelled tile
  fsm.TILE = TILE

  # Defaults different for forest tile
  if (fsm.TILE == "forest")
    fsm.hfsn = Tf(0.3)
    fsm.z0sn = Tf(0.01)
  end

  # Settings specific for DENSITY=0
  if (fsm.DENSTY == 0)
    fsm.rhof = fsm.rho0
  end

  # Initialize surface properties for non-default tiles
  if (fsm.TILE == "glacier")
    fsm.alb0[:, :] .= Tf(0.3)
    fsm.z0sf[:, :] .= Tf(0.04)
  end

  # Derived soil parameters
  for j = 1:fsm.Ny
    for i = 1:fsm.Nx
      if (fsm.fcly[i, j] + fsm.fsnd[i, j] > Tf(1))
        fsm.fcly[i, j] = Tf(1) - fsm.fsnd[i, j]
      end
      fsm.b[i, j] = Tf(3.1) + Tf(15.7) * fsm.fcly[i, j] - Tf(0.3) * fsm.fsnd[i, j]
      fsm.hcap_soil[i, j] = (Tf(2.128) * fsm.fcly[i, j] + Tf(2.385) * fsm.fsnd[i, j]) * Tf(1e6) / (fsm.fcly[i, j] + fsm.fsnd[i, j])
      fsm.sathh[i, j] = Tf(10)^(Tf(0.17) - Tf(0.63) * fsm.fcly[i, j] - Tf(1.58) * fsm.fsnd[i, j])
      fsm.Vsat[i, j] = Tf(0.505) - Tf(0.037) * fsm.fcly[i, j] - Tf(0.142) * fsm.fsnd[i, j]
      fsm.Vcrit[i, j] = fsm.Vsat[i, j] * (fsm.sathh[i, j] / Tf(3.364))^(Tf(1) / fsm.b[i, j])
      hcon_min = (hcon_clay^fsm.fcly[i, j]) * (hcon_sand^(Tf(1) - fsm.fcly[i, j]))
      fsm.hcon_soil[i, j] = (hcon_air^fsm.Vsat[i, j]) * (hcon_min^(Tf(1) - fsm.Vsat[i, j]))
    end
  end

  # Initial soil profiles
  fsat = Tf(0.5)
  Tprof = Tf(285)
  for k = 1:fsm.Nsoil
    fsm.theta[k, :, :] .= fsat * fsm.Vsat[:, :]
    fsm.Tsoil[k, :, :] .= Tprof
  end

  # Load terrain properties from landuse data
  fsm.fsky_terr .= Tf.(landuse["skyvf"]["data"])
  fsm.lat .= Tf.(landuse["y"]["data"])  # TODO remove
  fsm.lon .= Tf.(landuse["x"]["data"])  # TODO remove
  fsm.dem .= Tf.(landuse["dem"]["data"])

  # Cap surface temperatures for glacier 
  if (fsm.TILE == "glacier")
    for j = 1:fsm.Ny
      for i = 1:fsm.Nx
        fsm.Tsrf[i, j] = min(fsm.Tsrf[i, j], Tm)
        for k = 1:fsm.Nsoil
          fsm.Tsoil[k, i, j] = min(fsm.Tsoil[k, i, j], Tm)
        end
      end
    end
  end

  # Model tile fractions 
  if (fsm.TILE == "open")
    fsm.tilefrac = ones(Tf, size(fsm.dem))
  else
    fsm.tilefrac .= Tf.(landuse[lowercase(TILE)]["data"])
  end
  
  # Initialize snow cover fraction specific variables
  fsm.slopemu .= Tf.(landuse["slopemu"]["data"])
  fsm.xi .= Tf.(landuse["xi"]["data"])
  fsm.Ld .= Tf.(landuse["Ld"]["data"])

  # Canopy properties
  if (fsm.TILE != "forest")
    fsm.VAI[:, :] .= Tf(0)
    fsm.hcan[:, :] .= Tf(0)
    fsm.fsky[:, :] .= Tf(1)
    fsm.trcn[:, :] .= exp.(-fsm.kdif .* fsm.VAI[:, :])
    fsm.fveg[:, :] .= Tf(1) .- exp.(-fsm.kveg .* fsm.VAI[:, :])
    fsm.fves[:, :] .= Tf(1) .- exp.(-fsm.kveg .* fsm.VAI[:, :])
  else
    fsm.Qcan .= Tf(0)
    fsm.Sveg .= Tf(0)
    fsm.Tcan .= Tf(285)
    fsm.Tveg .= Tf(285)

    fsm.fveg .= Tf.(landuse["fveg"]["data"])
    fsm.hcan .= Tf.(landuse["hcan"]["data"])
    fsm.lai .= Tf.(landuse["lai"]["data"])
    fsm.vfhp .= Tf.(landuse["vfhp"]["data"])
    fsm.fves .= Tf.(landuse["fves"]["data"])
    fsm.pmultf .= Tf.(landuse["prec_multi"]["data"])

    fsm.VAI[:, :] = fsm.lai[:, :]
    fsm.trcn[:, :] = Tf(1) .- Tf(0.9) .* fsm.fveg[:, :]
    for j = 1:fsm.Ny
      for i = 1:fsm.Nx
        fsm.fsky[i, j] = fsm.vfhp[i, j] ./ fsm.trcn[i, j]
        if (fsm.fsky[i, j] > Tf(1))
          fsm.trcn[i, j] = fsm.vfhp[i, j]
        end
        if (fsm.fsky[i, j] > Tf(1))
          fsm.fsky[i, j] = Tf(1)
        end
      end
    end
  end

  fsm.canh[:, :] = Tf(12500) * fsm.VAI[:, :]
  fsm.scap[:, :] = fsm.cvai * fsm.VAI[:, :]

  # Tuned snow surface properties

  if fsm.OSHDTN == 0

    fsm.adm = Tf(100)
    fsm.adc[:, :] = Tf(1000)
    fsm.afs[:, :] = fsm.asmx
    if (fsm.TILE == "glacier" || ((fsm.SNTRAN == 1 || fsm.SNSLID == 1) && fsm.glacierfrac(i, j) > eps(Tf)))
      fsm.z0_snow[:, :] = Tf(0.0009)
    else
      fsm.z0_snow[:, :] = fsm.z0sn
    end

  else

    fsm.adm = Tf(130)

    for j = 1:fsm.Ny
      for i = 1:fsm.Nx

        # Elevation-dependent tuning of cold snow albedo decay time
        if (fsm.dem[i, j] >= Tf(2300))
          fsm.adc[i, j] = Tf(6000)
        elseif (fsm.dem[i, j] <= Tf(1500))
          fsm.adc[i, j] = Tf(3000)
        else
          fsm.adc[i, j] = Tf(6000) + (Tf(2300) - fsm.dem[i, j]) / (Tf(2300) - Tf(1500)) * (Tf(3000) - Tf(6000))
        end

        # Fresh snow albedo is now constant (previously elevation-dependent)
        fsm.afs[i, j] = fsm.asmx

        # Elevation-dependent tuning of snow roughness length
        if (fsm.TILE == "glacier" || ((fsm.SNTRAN == 1 || fsm.SNSLID == 1) && fsm.glacierfrac[i, j] > eps(Tf)))
          fsm.z0_snow[i, j] = Tf(0.0009)
        elseif (fsm.TILE == "forest")
          fsm.z0_snow[i, j] = fsm.z0sn
        else
          if (fsm.dem[i, j] >= Tf(2300))
            fsm.z0_snow[i, j] = Tf(0.01)
          elseif (fsm.dem[i, j] >= Tf(1500))
            fsm.z0_snow[i, j] = Tf(0.2) + (fsm.dem[i, j] - Tf(1500)) / (Tf(2300) - Tf(1500)) * (Tf(0.01) - Tf(0.2))
          else
            fsm.z0_snow[i, j] = Tf(0.2)
          end
        end

      end
    end

  end

  return fsm

end
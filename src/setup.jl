"""
    setup(Tf, Ti, landuse, Nx, Ny, settings)

Initialize the FSM snow model with specified configuration and domain properties.

Creates and configures the FSM model state structure with landuse data, model parameters,
and domain dimensions. Sets up soil properties, surface characteristics, and applies 
configuration-specific settings for different surface types and model behaviors.

# Arguments
- `Tf`: Floating-point precision type (typically Float32 or Float64)
- `Ti`: Integer type for array indices (typically Int32 or Int64)
- `landuse::Dict`: Landuse data dictionary with topographic and surface properties
- `Nx::Int, Ny::Int`: Model domain dimensions
- `settings::Dict`: Configuration dictionary containing:
  - `"tile"`: Surface tile type ("open", "forest", "glacier")
  - `"config"` (optional): Model configuration flags (SNFRAC, CANMOD, EXCHNG, ZOFFST, etc.)
  - `"params"` (optional): Parameter overrides (hfsn, z0sn, etc.)

# Returns
- `FSM`: Initialized model state structure ready for simulation
"""
function setup(Tf, Ti, landuse::Dict, Nx::Int, Ny::Int, settings::Dict)

  @unpack_constants(Tf)

  # Create fsm object
  fsm = FSM{Tf,Ti}(Nx=Nx, Ny=Ny)

  # Set tile
  fsm.TILE = settings["tile"]

  # Apply model configuration
  if haskey(settings, "config")
    for (key, value) in settings["config"]
      setfield!(fsm, Symbol(key), Int32(value))
    end
  end

  # Apply parameter overrides
  if haskey(settings, "params")
    for (key, value) in settings["params"]
      setfield!(fsm, Symbol(key), Tf(value))
    end
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

  mask = fsm.fcly .+ fsm.fsnd .> Tf(1)
  fsm.fcly[mask] .= Tf(1) .- fsm.fsnd[mask]
    
  fsm.b .= Tf(3.1) .+ Tf(15.7) .* fsm.fcly .- Tf(0.3) .* fsm.fsnd
  fsm.hcap_soil .= (Tf(2.128) .* fsm.fcly .+ Tf(2.385) .* fsm.fsnd) .* Tf(1e6) ./ (fsm.fcly .+ fsm.fsnd)
  fsm.sathh .= Tf(10) .^ (Tf(0.17) .- Tf(0.63) .* fsm.fcly .- Tf(1.58) .* fsm.fsnd)
  fsm.Vsat .= Tf(0.505) .- Tf(0.037) .* fsm.fcly .- Tf(0.142) .* fsm.fsnd
  fsm.Vcrit .= fsm.Vsat .* (fsm.sathh ./ Tf(3.364)) .^ (Tf(1) ./ fsm.b)
  hcon_min = (hcon_clay .^ fsm.fcly) .* (hcon_sand .^ (Tf(1) .- fsm.fcly))
  fsm.hcon_soil .= (hcon_air .^ fsm.Vsat) .* (hcon_min .^ (Tf(1) .- fsm.Vsat))

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
    fsm.Tsrf .= min.(fsm.Tsrf, Tm)
    fsm.Tsoil .= min.(fsm.Tsoil, Tm)
  end

  # Model tile fractions 
  if (fsm.TILE == "open")
    fsm.tilefrac = ones(Tf, size(fsm.dem))
  else
    fsm.tilefrac .= Tf.(landuse[lowercase(fsm.TILE)]["data"])
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
    fsm.fsky .= fsm.vfhp ./ fsm.trcn
    # Handle values where fsky > 1
    mask = fsm.fsky .> Tf(1)
    fsm.trcn[mask] .= fsm.vfhp[mask]
    fsm.fsky[mask] .= Tf(1)
  end

  fsm.canh[:, :] = Tf(12500) * fsm.VAI[:, :]
  fsm.scap[:, :] = fsm.cvai * fsm.VAI[:, :]

  # Tuned snow surface properties

  # if fsm.OSHDTN == 0

  #   fsm.adm = Tf(100)
  #   fsm.adc[:, :] .= Tf(1000)
  #   fsm.afs[:, :] .= fsm.asmx
  #   if (fsm.TILE == "glacier" || ((fsm.SNTRAN == 1 || fsm.SNSLID == 1) && fsm.glacierfrac(i, j) > eps(Tf)))
  #     fsm.z0_snow[:, :] .= Tf(0.0009)
  #   else
  #     fsm.z0_snow[:, :] .= fsm.z0sn
  #   end

  # else

    fsm.adm = Tf(130)

    # Elevation-dependent tuning of cold snow albedo decay time
    fsm.adc .= Tf(6000) .+ (Tf(2300) .- fsm.dem) ./ (Tf(2300) .- Tf(1500)) .* (Tf(3000) .- Tf(6000))
    fsm.adc[fsm.dem .>= Tf(2300)] .= Tf(6000)
    fsm.adc[fsm.dem .<= Tf(1500)] .= Tf(3000)

    for j = 1:fsm.Ny
      for i = 1:fsm.Nx

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

  # end

  fsm.z0_snow .= Tf(0.01)

  return fsm

end
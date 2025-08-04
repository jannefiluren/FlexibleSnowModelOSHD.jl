function setup_matfiles_point(Tf, Ti, landuse::Dict, Nx::Int, Ny::Int; SNFRAC = nothing)

  # Constants

  undefined = Tf(1.e+6)

  # # Read namelist file

  # f90nml = pyimport("f90nml")
  # nml = f90nml.read(joinpath(folder, "OPTIONS.nam"))
  # nml = pyconvert(Dict, nml)

  # Initialization of variables

  # nam_grid = pyconvert(Dict,nml["nam_grid"])

  Nsmax      = 3
  Nsoil      = 4
  # Nx         = 1   # TODO fix
  # Ny         = 1   # TODO fix
  Ds_min     = Tf(0.01)
  Ds_surflay = Tf(0.5)

  # Create fsm object

  fsm = FSM{Tf, Ti}(Nsmax=Nsmax, Nsoil=Nsoil, Nx=Nx, Ny=Ny)  # TODO add Ds_min and Ds_surflay

  # Driving data

  # nam_driving = pyconvert(Dict,nml["nam_driving"])

  fsm.zT   = Tf(10)
  fsm.zU   = Tf(10)
  fsm.zRH  = Tf(10)

  # Model configuration

  # nam_modconf = pyconvert(Dict,nml["nam_modconf"])

  fsm.ALBEDO = 2
  fsm.CANMOD = 0
  fsm.CONDCT = 1
  fsm.DENSTY = 3
  fsm.EXCHNG = 1
  fsm.HYDROL = 2
  fsm.SNFRAC = SNFRAC === nothing ? 3 : SNFRAC
  fsm.RADSBG = 0
  fsm.ZOFFST = 0
  fsm.OSHDTN = 1
  fsm.ALRADT = 0
  fsm.SNTRAN = 0
  fsm.SNSLID = 0
  fsm.SNOLAY = 0
  fsm.CHECKS = 0
  fsm.HN_ON  = false
  fsm.FOR_HN = true

  # Model perturbations

  # nam_modpert = pyconvert(Dict,nml["nam_modpert"])

  fsm.Z0PERT = false
  fsm.WCPERT = false
  fsm.FSPERT = false
  fsm.ALPERT = false
  fsm.SLPERT = false

  # Modelled tile

  # nam_modtile = pyconvert(Dict,nml["nam_modtile"])

  fsm.TILE = "open"
  fsm.tthresh = Tf(0.1)

  # Defaults for numerical solution parameters
  fsm.Nitr = 4

  # Defaults for canopy parameters
  fsm.avg0 = Tf(0.1)
  fsm.avgs = Tf(0.4)
  fsm.cden = Tf(0.004)
  fsm.cvai = Tf(4.4)
  fsm.cveg = Tf(20)
  fsm.Gcn1 = Tf(0.5)
  fsm.Gcn2 = Tf(0)
  fsm.gsnf = Tf(0)
  fsm.kdif = Tf(0.5)
  fsm.kveg = Tf(1)
  fsm.rchd = Tf(0.67)
  fsm.rchz = Tf(0.2) 
  fsm.tcnc = Tf(240)
  fsm.tcnm = Tf(48)

  # Defaults for snow parameters
  fsm.a_eta = Tf(0.1)
  fsm.asmx = Tf(0.86)       # unused if OSHDTN = 1
  fsm.asmn = Tf(0.6)
  fsm.b_eta = Tf(0.023)
  fsm.bstb = Tf(5)
  fsm.bthr = Tf(2)
  fsm.c_eta = Tf(250)
  fsm.eta0 = Tf(3.7e7)
  fsm.eta1 = Tf(7.62237e6)
  fsm.hfsn = Tf(0.1)
  fsm.kfix = Tf(0.24)
  fsm.rho0 = Tf(300)
  fsm.rhob = Tf(6)
  fsm.rhoc = Tf(26)
  fsm.rhof = Tf(109)
  fsm.rhos_min = Tf(50.0)
  fsm.rhos_max = Tf(750.0)
  fsm.rcld = Tf(300)
  fsm.rgr0 = Tf(5e-5)
  fsm.rmlt = Tf(500)
  fsm.snda = Tf(2.8e-6)
  fsm.Talb = Tf(-2)
  fsm.tcld = Tf(1000)
  fsm.tmlt = Tf(100)
  fsm.trho = Tf(200)
  fsm.Wirr = Tf(0.03)
  fsm.z0sn = Tf(0.002)
  fsm.Sfmin = Tf(10)

  # some defaults different for forest tile - commented-out values based on FS-EBP runs, revisit during tuning
  if (fsm.TILE == "forest")
    # asmx = Tf(0.88)
    fsm.hfsn = Tf(0.3)    
    fsm.z0sn = Tf(0.01)   
  end

  # Defaults for ground surface parameters
  fsm.bstb = Tf(5)
  fsm.gsat = Tf(0.01)

  # Defaults for additional parameters required for forest snow process parametrization
  fsm.adfs = Tf(3)
  fsm.adfl = Tf(2)
  fsm.fsar = Tf(0.1)
  fsm.psf  = Tf(1)
  fsm.psr  = Tf(0.1)
  fsm.wcan = Tf(2.5)
  fsm.zsub = Tf(2)
  fsm.zgf = Tf(1)
  fsm.zgr = Tf(0)
  fsm.khcf = Tf(3)

  if (fsm.DENSTY == 0)
    fsm.rhof = fsm.rho0
  end

  # Surface data from defaults or namelists
  # Surface properties 
  # if (SNTRAN == 1) allocate(vegsnowd_xy(Nx,Ny))    TODO: translate later
  if (fsm.TILE == "glacier")
    fsm.alb0[:,:] .= Tf(0.3)
    fsm.z0sf[:,:] .= Tf(0.04)
  else
    fsm.alb0[:,:] .= Tf(0.2)
    fsm.z0sf[:,:] .= Tf(0.2)
  end
  fsm.fcly[:,:] .= Tf(0.3)
  fsm.fsnd[:,:] .= Tf(0.6)
  if (fsm.TILE == "forest")
    fsm.z0sf[:,:] .= Tf(0.2)
  end
  if (fsm.SNTRAN == 1)
    fsm.vegsnowd_xy[:,:] .= Tf(0.1)
  end

  # Canopy parameters
  fsm.canh[:,:]        .= undefined
  fsm.fsky[:,:]        .= undefined
  fsm.fveg[:,:]        .= undefined
  fsm.fves[:,:]        .= undefined
  fsm.hcan[:,:]        .= undefined
  fsm.pmultf[:,:]      .= undefined
  fsm.scap[:,:]        .= undefined
  fsm.trcn[:,:]        .= undefined
  fsm.VAI[:,:]         .= undefined

  # Terrain properties
  fsm.fsky_terr[:,:]   .= undefined
  fsm.lat[:,:]         .= undefined
  fsm.lon[:,:]         .= undefined
  fsm.dem[:,:]         .= undefined
  fsm.tilefrac[:,:]    .= undefined
  fsm.glacierfrac[:,:] .= undefined

  if (fsm.SNFRAC == 0)
    fsm.slopemu[:,:]   .= undefined
    fsm.xi[:,:]        .= undefined
  end
  if (fsm.SNFRAC == 0 || fsm.SNTRAN == 1)
    fsm.Ld[:,:]        .= undefined
  end

  if (fsm.SNSLID == 1)
    fsm.slope[:,:]     .= undefined
    fsm.Shd[:,:]       .= undefined
  end

  # Derived soil parameters
  for j = 1:fsm.Ny
    for i = 1:fsm.Nx
      if (fsm.fcly[i,j] + fsm.fsnd[i,j] > Tf(1))
        fsm.fcly[i,j] = Tf(1) - fsm.fsnd[i,j]
      end
      fsm.b[i,j] = Tf(3.1) + Tf(15.7)*fsm.fcly[i,j] - Tf(0.3)*fsm.fsnd[i,j]
      fsm.hcap_soil[i,j] = (Tf(2.128)*fsm.fcly[i,j] + Tf(2.385)*fsm.fsnd[i,j])*Tf(1e6) / (fsm.fcly[i,j] + fsm.fsnd[i,j])
      fsm.sathh[i,j] = Tf(10)^(Tf(0.17) - Tf(0.63)*fsm.fcly[i,j] - Tf(1.58)*fsm.fsnd[i,j])
      fsm.Vsat[i,j] = Tf(0.505) - Tf(0.037)*fsm.fcly[i,j] - Tf(0.142)*fsm.fsnd[i,j]
      fsm.Vcrit[i,j] = fsm.Vsat[i,j]*(fsm.sathh[i,j]/Tf(3.364))^(Tf(1)/fsm.b[i,j])
      hcon_min = (hcon_clay^fsm.fcly[i,j]) * (hcon_sand^(Tf(1) - fsm.fcly[i,j]))
      fsm.hcon_soil[i,j] = (hcon_air^fsm.Vsat[i,j]) * (hcon_min^(Tf(1) - fsm.Vsat[i,j]))
    end
  end

  # Convert time scales from hours to seconds
  ### fsm.dt = 3600*fsm.dt    # TODO this is already done in types.jl - check where to to this transformation in future
  fsm.tcnc = Tf(3600)*fsm.tcnc
  fsm.tcnm = Tf(3600)*fsm.tcnm
  fsm.tcld = Tf(3600)*fsm.tcld
  fsm.tmlt = Tf(3600)*fsm.tmlt
  fsm.trho = Tf(3600)*fsm.trho

  # Default initialization of state variables 
  fsm.albs[:,:]               .= undefined
  fsm.Ds[:,:,:]               .= undefined
  fsm.fsnow[:,:]              .= undefined
  fsm.Nsnow[:,:]              .= undefined
  fsm.Qcan[:,:]               .= undefined
  fsm.rgrn[:,:,:]             .= undefined #*GM watch out: rgrn currently not tracked
  fsm.histowet[:,:,:]         .= undefined # LQ: histowet only tracked if SNTRAN
  fsm.Sice[:,:,:]             .= undefined
  fsm.Sliq[:,:,:]             .= undefined
  fsm.Sveg[:,:]               .= undefined
  fsm.Tcan[:,:]               .= undefined
  fsm.Tsnow[:,:,:]            .= undefined
  fsm.Tsoil[:,:,:]            .= undefined
  fsm.Tveg[:,:]               .= undefined

  # Allocation and initialization of optional state variables
  if (fsm.SNFRAC == 0 || fsm.SNFRAC == 2)
    fsm.snowdepthmax[:,:] .= undefined
  end

  if (fsm.SNFRAC == 0)
    fsm.snowdepthmin[:,:]     .= undefined
    fsm.snowdepthhist[:,:,:]  .= undefined
    fsm.swemin[:,:]           .= undefined
    fsm.swemax[:,:]           .= undefined
    fsm.swehist[:,:,:]        .= undefined
  end

  if (fsm.SNTRAN == 1)
    fsm.dSWE_tot_subl[:,:]    .= undefined
    fsm.dSWE_tot_salt[:,:]    .= undefined
    fsm.dSWE_tot_susp[:,:]    .= undefined
  end

  if (fsm.SNSLID == 1)
    dSWE_tot_slide[:,:]   .= undefined
    index_sorted_dem[:,:] .= undefined
  end

  # Initial soil profiles from namelist
  fsat  = Tf(0.5)
  Tprof = Tf(285)
  for k = 1:fsm.Nsoil
    fsm.theta[k,:,:] .= fsat*fsm.Vsat[:,:]
    fsm.Tsoil[k,:,:] .= Tprof
  end
  fsm.Tsrf[:,:] = fsm.Tsoil[1,:,:]


  # Code inserted from OPEN_FILES.F90

  # L220-224
  if (fsm.TILE == "forest")
    landuse_tile             = "landuse_forest.bin"
  else (fsm.TILE == "glacier")
    landuse_tile             = "landuse_glacier.bin"
  end

  #L 225-227
  if (fsm.SNTRAN == 1 || fsm.SNSLID == 1)
    landuse_glac             = "landuse_glacier.bin"
  end


  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #
  ##-4- !!!!!!!!!!!!!!!!!!!! READ DRIVING/STATES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #! states relevant to both open and forest simulation
  fsm.albs .= Tf(0.85)                                                               # read!(joinpath(folder, "states_in_alse.bin"), fsm.albs)
  fsm.Ds .= Tf(0)                                                                    # read!(joinpath(folder, "states_in_hsnl.bin"), fsm.Ds)
  fsm.fsnow .= Tf(0)                                                                 # read!(joinpath(folder, "states_in_scfe.bin"), fsm.fsnow)
  fsm.Nsnow .= Tf(0)                                                                 # read!(joinpath(folder, "states_in_nsne.bin"), fsm.Nsnow)
  fsm.Sice .= Tf(0)                                                                  # read!(joinpath(folder, "states_in_sicl.bin"), fsm.Sice)
  fsm.Sliq .= Tf(0)                                                                  # read!(joinpath(folder, "states_in_slql.bin"), fsm.Sliq)
  fsm.Tsrf .= 285                                                                     # read!(joinpath(folder, "states_in_tsfe.bin"), fsm.Tsrf)
  fsm.Tsnow .= Tm                                                                    # read!(joinpath(folder, "states_in_tsnl.bin"), fsm.Tsnow)
  fsm.Tsoil .= Tf(285)                                                               # read!(joinpath(folder, "states_in_tsll.bin"), fsm.Tsoil)
  fsm.fsky_terr .= Tf.(landuse["skyvf"]["data"][landuse["is_domain"]])               # read!(joinpath(folder, "landuse_skyvf.bin"), fsm.fsky_terr)
  fsm.lat .= Tf.(landuse["y"][landuse["is_domain"]])                                 # read!(joinpath(folder, "landuse_lat.bin"), fsm.lat)   # TODO this wrong - currently swiss coords
  fsm.lon .= Tf.(landuse["x"][landuse["is_domain"]])                                 # read!(joinpath(folder, "landuse_lon.bin"), fsm.lon)   # TODO this wrong - currently swiss coords
  fsm.dem .= Tf.(landuse["dem"]["data"][landuse["is_domain"]])                       # read!(joinpath(folder, "landuse_dem.bin"), fsm.dem)

  # Cap glacier temperatures to 0°C 
  if (fsm.TILE == "glacier")
    for j = 1:fsm.Ny
    for i = 1:fsm.Nx
      fsm.Tsrf[i,j] = min(fsm.Tsrf[i,j],Tm)
      for k = 1:fsm.Nsoil
        fsm.Tsoil[k,i,j] = min(fsm.Tsoil[k,i,j],Tm)
      end
    end
    end
  end

  # model tile fractions 
  if (fsm.TILE == "open")
    fsm.tilefrac = fsm.dem./fsm.dem   # temporary fix to get ones within our entire domain, assuming we always want to run an open tile. may have to be revisited
    if (fsm.SNTRAN == 1 || fsm.SNSLID == 1)
      # read!(joinpath(folder, landuse_glac), fsm.glacierfrac)    TODO add this later
    end
  else
    # read!(joinpath(folder, landuse_tile), fsm.tilefrac)   TODO add this later
  end 

  if (fsm.TILE == "open" && (fsm.SNTRAN == 1 || fsm.SNSLID == 1))   # Cap glacier temperatures only in relevant pixels, and alter snow free albedo and roughness lengths
    for j = 1:Ny
    for i = 1:Nx
      if (fsm.glacierfrac[i,j] > eps(Tf))
        fsm.Tsrf[i,j] = min(fsm.Tsrf[i,j],Tm)
        fsm.alb0[i,j] = Tf(0.3) 
        fsm.z0sf[i,j] = Tf(0.04)
        for k = 1:fsm.Nsoil
          fsm.Tsoil[k,i,j] = min(fsm.Tsoil[k,i,j],Tm)
        end
      end
    end
    end
  end

  if (fsm.SNFRAC == 0 || fsm.SNFRAC == 2)
    fsm.snowdepthmax .= 0     # read!(joinpath(folder, "states_in_hsmx.bin"), fsm.snowdepthmax)
  end

  if (fsm.SNFRAC == 0)
    # states specific to open runs
    fsm.snowdepthmin .= Tf(0)  # read!(joinpath(folder, "states_in_hsmn.bin"), fsm.snowdepthmin)
    fsm.snowdepthhist .= Tf(0) # read!(joinpath(folder, "states_in_hshs.bin"), fsm.snowdepthhist)
    fsm.swemin .= Tf(0)   # read!(joinpath(folder, "states_in_swmn.bin"), fsm.swemin)
    fsm.swemax .= Tf(0)   # read!(joinpath(folder, "states_in_swmx.bin"), fsm.swemax)
    fsm.swehist .= Tf(0)   # read!(joinpath(folder, "states_in_swhs.bin"), fsm.swehist)
    fsm.slopemu .= Tf.(landuse["slopemu"][landuse["is_domain"]])   # read!(joinpath(folder, "landuse_slopemu.bin"), fsm.slopemu)
    fsm.xi .= Tf.(landuse["xi"][landuse["is_domain"]])      # read!(joinpath(folder, "landuse_xi.bin"), fsm.xi)
  end

  if (fsm.SNFRAC == 0 || fsm.SNTRAN == 1)
    fsm.Ld .= Tf.(landuse["Ld"][landuse["is_domain"]])    # read!(joinpath(folder, "landuse_Ld.bin"), fsm.Ld)
  end

  if (fsm.TILE != "forest")
    # canopy properties (no canopy)
    fsm.VAI[:,:]  .= Tf(0)
    fsm.hcan[:,:] .= Tf(0)
    fsm.fsky[:,:] .= Tf(1)
    fsm.trcn[:,:] .= exp.(-fsm.kdif.*fsm.VAI[:,:])
    fsm.fveg[:,:] .= Tf(1) .- exp.(-fsm.kveg.*fsm.VAI[:,:])
    fsm.fves[:,:] .= Tf(1) .- exp.(-fsm.kveg.*fsm.VAI[:,:])
  else # TILE == 'forest'
    # lus fields specific to forest runs
    # read!(joinpath(folder, "states_in_qcan.bin"), fsm.Qcan)  TODO add later
    # read!(joinpath(folder, "states_in_sveg.bin"), fsm.Sveg)
    # read!(joinpath(folder, "states_in_tcan.bin"), fsm.Tcan)
    # read!(joinpath(folder, "states_in_tveg.bin"), fsm.Tveg)
    # read!(joinpath(folder, "landuse_fveg.bin"), fsm.fveg)
    # read!(joinpath(folder, "landuse_hcan.bin"), fsm.hcan)
    # read!(joinpath(folder, "landuse_lai.bin"), fsm.lai)
    # read!(joinpath(folder, "landuse_vfhp.bin"), fsm.vfhp)
    # read!(joinpath(folder, "landuse_fves.bin"), fsm.fves)
    # read!(joinpath(folder, "landuse_pmultf.bin"), fsm.pmultf)
    
    # derived canopy properties 
    fsm.VAI[:,:] = fsm.lai[:,:]
    trcn[:,:] = Tf(1) .- Tf(0.9).*fsm.fveg[:,:]  
    for j = 1:fsm.Ny
      for i = 1:fsm.Nx
        fsm.fsky[i,j] = fsm.vfhp[i,j]./fsm.trcn[i,j]
        if ( fsm.fsky[i,j] > Tf(1) )
          fsm.trcn[i,j] = fsm.vfhp[i,j]
        end
        if ( fsm.fsky[i,j] > Tf(1) )
          fsm.fsky[i,j] = Tf(1)
        end
      end
    end 
  end

  # derived canopy parameters
  fsm.canh[:,:] = Tf(12500)*fsm.VAI[:,:]
  fsm.scap[:,:] = fsm.cvai*fsm.VAI[:,:]

  if (fsm.SNTRAN == 1)
    # states specific to SNOWTRAN3D
    # read!(joinpath(folder, "states_in_hiwl.bin"), histowet)  TODO add later
    # read!(joinpath(folder, "states_in_sblt.bin"), dSWE_tot_subl)
    # read!(joinpath(folder, "states_in_sltt.bin"), dSWE_tot_salt)
    # read!(joinpath(folder, "states_in_sspt.bin"), dSWE_tot_susp)
  end

  if (fsm.SNSLID == 1)
    # states specific to SnowSlide
    # read!(joinpath(folder, "landuse_slope.bin"), slope)  TODO add later
    # read!(joinpath(folder, "landuse_shd.bin"), Shd)
    # read!(joinpath(folder, "states_in_sldt.bin"), dSWE_tot_slide)
    # read!(joinpath(folder, "states_in_idem.bin"), index_sorted_dem)
  end

  # Tuned snow surface properties

  if fsm.OSHDTN == 0

    fsm.adm = Tf(100)
    fsm.adc[:,:] = Tf(1000)
    fsm.afs[:,:] = fsm.asmx
    if (fsm.TILE == "glacier" || ((fsm.SNTRAN == 1 || fsm.SNSLID == 1) && fsm.glacierfrac(i,j) > eps(Tf)))
      fsm.z0_snow[:,:] = Tf(0.0009)
    else
      fsm.z0_snow[:,:] = fsm.z0sn
    end

  else

    fsm.adm = Tf(130)

    for j = 1:fsm.Ny
      for i = 1:fsm.Nx
      
        # Elevation-dependent tuning of cold snow albedo decay time
        if (fsm.dem[i,j] >= Tf(2300))
          fsm.adc[i,j]  = Tf(6000)
        elseif (fsm.dem[i,j] <= Tf(1500))
          fsm.adc[i,j] = Tf(3000)
        else
          fsm.adc[i,j] = Tf(6000) + (Tf(2300) - fsm.dem[i,j]) / (Tf(2300) - Tf(1500)) * (Tf(3000) - Tf(6000))
        end
        
        # Fresh snow albedo is now constant (previously elevation-dependent)
        fsm.afs[i,j] = fsm.asmx
        
        # Elevation-dependent tuning of snow roughness length
        if (fsm.TILE == "glacier" || ((fsm.SNTRAN == 1 || fsm.SNSLID == 1) && fsm.glacierfrac[i,j] > eps(Tf)))
          fsm.z0_snow[i,j] = Tf(0.0009)
        elseif (fsm.TILE == "forest")
          fsm.z0_snow[i,j] = fsm.z0sn
        else
          if (fsm.dem[i,j] >= Tf(2300))
            fsm.z0_snow[i,j] = Tf(0.01)
          elseif (fsm.dem[i,j] >= Tf(1500))
            fsm.z0_snow[i,j] = Tf(0.2) + (fsm.dem[i,j] - Tf(1500)) / (Tf(2300) - Tf(1500)) * (Tf(0.01) - Tf(0.2))
          else
            fsm.z0_snow[i,j] = Tf(0.2)
          end
        end
        
      end
    end

  end

  return fsm

end
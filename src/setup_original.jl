function setup_original(Tf, Ti, folder::String)

  # Constants

  undefined = 1.e+6

  # Read namelist file

  f90nml = pyimport("f90nml")
  nml = f90nml.read(joinpath(folder, "OPTIONS.nam"))
  nml = pyconvert(Dict, nml)

  # Initialization of variables

  nam_grid = pyconvert(Dict,nml["nam_grid"])

  Nsmax      = haskey(nam_grid, lowercase("NNsmax")) ? nam_grid[lowercase("NNsmax")] : 3
  Nsoil      = haskey(nam_grid, lowercase("NNsoil")) ? nam_grid[lowercase("NNsoil")] : 4
  Nx         = haskey(nam_grid, lowercase("NNx")) ? nam_grid[lowercase("NNx")] : 1
  Ny         = haskey(nam_grid, lowercase("NNy")) ? nam_grid[lowercase("NNy")] : 1
  Ds_min     = haskey(nam_grid, lowercase("DDs_min")) ? nam_grid[lowercase("DDs_min")] : 0.01
  Ds_surflay = haskey(nam_grid, lowercase("DDs_surflay")) ? nam_grid[lowercase("DDs_surflay")] : 0.5

  # Create fsm object

  fsm = FSM{Tf, Ti}(Nsmax=Nsmax, Nsoil=Nsoil, Nx=Nx, Ny=Ny)  # TODO add Ds_min and Ds_surflay

  # Driving data

  nam_driving = pyconvert(Dict,nml["nam_driving"])

  fsm.zT     = haskey(nam_driving, lowercase("zzT")) ? nam_driving[lowercase("zzT")] : 10
  fsm.zU     = haskey(nam_driving, lowercase("zzU")) ? nam_driving[lowercase("zzU")] : 10
  fsm.zRH    = haskey(nam_driving, lowercase("zzRH")) ? nam_driving[lowercase("zzRH")] : 10

  # Model configuration

  nam_modconf = pyconvert(Dict,nml["nam_modconf"])

  fsm.ALBEDO = haskey(nam_modconf, lowercase("NALBEDO")) ? nam_modconf[lowercase("NALBEDO")] : -1
  fsm.CANMOD = haskey(nam_modconf, lowercase("NCANMOD")) ? nam_modconf[lowercase("NCANMOD")] : -1
  fsm.CONDCT = haskey(nam_modconf, lowercase("NCONDCT")) ? nam_modconf[lowercase("NCONDCT")] : -1
  fsm.DENSTY = haskey(nam_modconf, lowercase("NDENSTY")) ? nam_modconf[lowercase("NDENSTY")] : -1
  fsm.EXCHNG = haskey(nam_modconf, lowercase("NEXCHNG")) ? nam_modconf[lowercase("NEXCHNG")] : -1
  fsm.HYDROL = haskey(nam_modconf, lowercase("NHYDROL")) ? nam_modconf[lowercase("NHYDROL")] : -1
  fsm.SNFRAC = haskey(nam_modconf, lowercase("NSNFRAC")) ? nam_modconf[lowercase("NSNFRAC")] : -1
  fsm.RADSBG = haskey(nam_modconf, lowercase("NRADSBG")) ? nam_modconf[lowercase("NRADSBG")] : -1
  fsm.ZOFFST = haskey(nam_modconf, lowercase("NZOFFST")) ? nam_modconf[lowercase("NZOFFST")] : -1
  fsm.OSHDTN = haskey(nam_modconf, lowercase("NOSHDTN")) ? nam_modconf[lowercase("NOSHDTN")] : -1
  fsm.ALRADT = haskey(nam_modconf, lowercase("NALRADT")) ? nam_modconf[lowercase("NALRADT")] : -1
  fsm.SNTRAN = haskey(nam_modconf, lowercase("NSNTRAN")) ? nam_modconf[lowercase("NSNTRAN")] : -1
  fsm.SNSLID = haskey(nam_modconf, lowercase("NSNSLID")) ? nam_modconf[lowercase("NSNSLID")] : -1
  fsm.SNOLAY = haskey(nam_modconf, lowercase("NSNOLAY")) ? nam_modconf[lowercase("NSNOLAY")] : -1
  fsm.CHECKS = haskey(nam_modconf, lowercase("NCHECKS")) ? nam_modconf[lowercase("NCHECKS")] : -1
  fsm.HN_ON = haskey(nam_modconf, lowercase("LHN_ON")) ? nam_modconf[lowercase("LHN_ON")] : false
  fsm.FOR_HN = haskey(nam_modconf, lowercase("LFOR_HN")) ? nam_modconf[lowercase("LFOR_HN")] : false

  # Model perturbations

  nam_modpert = pyconvert(Dict,nml["nam_modpert"])

  fsm.Z0PERT = haskey(nam_modpert, lowercase("NZ0PERT")) ? nam_modpert[lowercase("NZ0PERT")] : false
  fsm.WCPERT = haskey(nam_modpert, lowercase("NWCPERT")) ? nam_modpert[lowercase("NWCPERT")] : false
  fsm.FSPERT = haskey(nam_modpert, lowercase("NFSPERT")) ? nam_modpert[lowercase("NFSPERT")] : false
  fsm.ALPERT = haskey(nam_modpert, lowercase("NALPERT")) ? nam_modpert[lowercase("NALPERT")] : false
  fsm.SLPERT = haskey(nam_modpert, lowercase("NSLPERT")) ? nam_modpert[lowercase("NSLPERT")] : false

  # Modelled tile

  nam_modtile = pyconvert(Dict,nml["nam_modtile"])

  fsm.TILE = haskey(nam_modtile, lowercase("CTILE")) ? nam_modtile[lowercase("CTILE")] : "open"
  fsm.tthresh = haskey(nam_modtile, lowercase("rtthresh")) ? nam_modtile[lowercase("rtthresh")] : 0.1

  # Defaults for numerical solution parameters
  fsm.Nitr = 4

  # Defaults for canopy parameters
  fsm.avg0 = 0.1
  fsm.avgs = 0.4
  fsm.cden = 0.004
  fsm.cvai = 4.4
  fsm.cveg = 20
  fsm.Gcn1 = 0.5
  fsm.Gcn2 = 0
  fsm.gsnf = 0
  fsm.kdif = 0.5
  fsm.kveg = 1
  fsm.rchd = 0.67
  fsm.rchz = 0.2 
  fsm.tcnc = 240
  fsm.tcnm = 48

  # Defaults for snow parameters
  fsm.a_eta = 0.1
  fsm.asmx = 0.86       # unused if OSHDTN = 1
  fsm.asmn = 0.6
  fsm.b_eta = 0.023
  fsm.bstb = 5
  fsm.bthr = 2
  fsm.c_eta = 250
  fsm.eta0 = 3.7e7
  fsm.eta1 = 7.62237e6
  fsm.hfsn = 0.1
  fsm.kfix = 0.24
  fsm.rho0 = 300
  fsm.rhob = 6
  fsm.rhoc = 26
  fsm.rhof = 109
  fsm.rhos_min = 50.0
  fsm.rhos_max = 750.0
  fsm.rcld = 300
  fsm.rgr0 = 5e-5
  fsm.rmlt = 500
  fsm.snda = 2.8e-6
  fsm.Talb = -2
  fsm.tcld = 1000
  fsm.tmlt = 100
  fsm.trho = 200
  fsm.Wirr = 0.03
  fsm.z0sn = 0.002
  fsm.Sfmin = 10

  # some defaults different for forest tile - commented-out values based on FS-EBP runs, revisit during tuning
  if (fsm.TILE == "forest")
    # asmx = 0.88
    fsm.hfsn = 0.3    
    fsm.z0sn = 0.01   
  end

  # Defaults for ground surface parameters
  fsm.bstb = 5
  fsm.gsat = 0.01

  # Defaults for additional parameters required for forest snow process parametrization
  fsm.adfs = 3
  fsm.adfl = 2
  fsm.fsar = 0.1
  fsm.psf  = 1
  fsm.psr  = 0.1
  fsm.wcan = 2.5
  fsm.zsub = 2
  fsm.zgf = 1
  fsm.zgr = 0
  fsm.khcf = 3

  if (fsm.DENSTY == 0)
    fsm.rhof = fsm.rho0
  end

  # Surface data from defaults or namelists
  # Surface properties 
  # if (SNTRAN == 1) allocate(vegsnowd_xy(Nx,Ny))    TODO: translate later
  if (fsm.TILE == "glacier")
    fsm.alb0[:,:] .= 0.3
    fsm.z0sf[:,:] .= 0.04
  else
    fsm.alb0[:,:] .= 0.2
    fsm.z0sf[:,:] .= 0.2
  end
  fsm.fcly[:,:] .= 0.3
  fsm.fsnd[:,:] .= 0.6
  if (fsm.TILE == "forest")
    fsm.z0sf[:,:] .= 0.2
  end
  if (fsm.SNTRAN == 1)
    fsm.vegsnowd_xy[:,:] .= 0.1
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
  fsm.tcnc = 3600*fsm.tcnc
  fsm.tcnm = 3600*fsm.tcnm
  fsm.tcld = 3600*fsm.tcld
  fsm.tmlt = 3600*fsm.tmlt
  fsm.trho = 3600*fsm.trho

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
    snowdepthmin[:,:]     .= undefined
    snowdepthhist[:,:,:]  .= undefined
    swemin[:,:]           .= undefined
    swemax[:,:]           .= undefined
    swehist[:,:,:]        .= undefined
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
  fsat  = 0.5
  Tprof = 285
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
  read!(joinpath(folder, "states_in_alse.bin"), fsm.albs)
  read!(joinpath(folder, "states_in_hsnl.bin"), fsm.Ds)
  read!(joinpath(folder, "states_in_scfe.bin"), fsm.fsnow)
  read!(joinpath(folder, "states_in_nsne.bin"), fsm.Nsnow)
  read!(joinpath(folder, "states_in_sicl.bin"), fsm.Sice)
  read!(joinpath(folder, "states_in_slql.bin"), fsm.Sliq)
  read!(joinpath(folder, "states_in_tsfe.bin"), fsm.Tsrf)
  read!(joinpath(folder, "states_in_tsnl.bin"), fsm.Tsnow)
  read!(joinpath(folder, "states_in_tsll.bin"), fsm.Tsoil)
  read!(joinpath(folder, "landuse_skyvf.bin"), fsm.fsky_terr)
  read!(joinpath(folder, "landuse_lat.bin"), fsm.lat)
  read!(joinpath(folder, "landuse_lon.bin"), fsm.lon)
  read!(joinpath(folder, "landuse_dem.bin"), fsm.dem)

  # Cap glacier temperatures to 0°C 
  if (fsm.TILE == "glacier")
    for j = 1:fsm.Ny
    for i = 1:fsm.Nx
      fsm.Tsrf[i,j] = min(fsm.Tsrf[i,j],Tm)
      for k = 1:fsm.Nsoil
        Tsoil[k,i,j] = min(Tsoil[k,i,j],Tm)
      end
    end
    end
  end

  # model tile fractions 
  if (fsm.TILE == "open")
    fsm.tilefrac .= fsm.dem./fsm.dem   # temporary fix to get ones within our entire domain, assuming we always want to run an open tile. may have to be revisited
    if (fsm.SNTRAN == 1 || fsm.SNSLID == 1)
      read!(joinpath(folder, landuse_glac), fsm.glacierfrac)
    end
  else
    read!(joinpath(folder, landuse_tile), fsm.tilefrac)
  end 

  if (fsm.TILE == "open" && (fsm.SNTRAN == 1 || fsm.SNSLID == 1))   # Cap glacier temperatures only in relevant pixels, and alter snow free albedo and roughness lengths
    for j = 1:Ny
    for i = 1:Nx
      if (fsm.glacierfrac[i,j] > eps(Tf))
        fsm.Tsrf[i,j] = min(fsm.Tsrf[i,j],Tm)
        fsm.alb0[i,j] = 0.3 
        fsm.z0sf[i,j] = 0.04
        for k = 1:fsm.Nsoil
          fsm.Tsoil[k,i,j] = min(fsm.Tsoil[k,i,j],Tm)
        end
      end
    end
    end
  end

  if (fsm.SNFRAC == 0 || fsm.SNFRAC == 2)
    read!(joinpath(folder, "states_in_hsmx.bin"), fsm.snowdepthmax)
  end

  if (fsm.SNFRAC == 0)
    # states specific to open runs
    read!(joinpath(folder, "states_in_hsmn.bin"), fsm.snowdepthmin)
    read!(joinpath(folder, "states_in_hshs.bin"), fsm.snowdepthhist)
    read!(joinpath(folder, "states_in_swmn.bin"), fsm.swemin)
    read!(joinpath(folder, "states_in_swmx.bin"), fsm.swemax)
    read!(joinpath(folder, "states_in_swhs.bin"), fsm.swehist)
    read!(joinpath(folder, "landuse_slopemu.bin"), fsm.slopemu)
    read!(joinpath(folder, "landuse_xi.bin"), fsm.xi)
  end

  if (fsm.SNFRAC == 0 || fsm.SNTRAN == 1)
    read!(joinpath(folder, "landuse_Ld.bin"), fsm.Ld)
  end

  if (fsm.TILE != "forest")
    # canopy properties (no canopy)
    fsm.VAI[:,:]  .= 0
    fsm.hcan[:,:] .= 0
    fsm.fsky[:,:] .= 1
    fsm.trcn[:,:] .= exp.(-fsm.kdif.*fsm.VAI[:,:])
    fsm.fveg[:,:] .= 1 .- exp.(-fsm.kveg.*fsm.VAI[:,:])
    fsm.fves[:,:] .= 1 .- exp.(-fsm.kveg.*fsm.VAI[:,:])
  else # TILE == 'forest'
    # lus fields specific to forest runs
    read!(joinpath(folder, "states_in_qcan.bin"), fsm.Qcan)
    read!(joinpath(folder, "states_in_sveg.bin"), fsm.Sveg)
    read!(joinpath(folder, "states_in_tcan.bin"), fsm.Tcan)
    read!(joinpath(folder, "states_in_tveg.bin"), fsm.Tveg)
    read!(joinpath(folder, "landuse_fveg.bin"), fsm.fveg)
    read!(joinpath(folder, "landuse_hcan.bin"), fsm.hcan)
    read!(joinpath(folder, "landuse_lai.bin"), fsm.lai)
    read!(joinpath(folder, "landuse_vfhp.bin"), fsm.vfhp)
    read!(joinpath(folder, "landuse_fves.bin"), fsm.fves)
    read!(joinpath(folder, "landuse_pmultf.bin"), fsm.pmultf)
    
    # derived canopy properties 
    fsm.VAI[:,:] = fsm.lai[:,:]
    trcn[:,:] = 1 .- 0.9.*fsm.fveg[:,:]  
    for j = 1:fsm.Ny
      for i = 1:fsm.Nx
        fsm.fsky[i,j] = fsm.vfhp[i,j]./fsm.trcn[i,j]
        if ( fsm.fsky[i,j] > 1 )
          fsm.trcn[i,j] = fsm.vfhp[i,j]
        end
        if ( fsm.fsky[i,j] > 1 )
          fsm.fsky[i,j] = 1
        end
      end
    end 
  end

  # derived canopy parameters
  fsm.canh[:,:] = 12500*fsm.VAI[:,:]
  fsm.scap[:,:] = fsm.cvai*fsm.VAI[:,:]

  if (fsm.SNTRAN == 1)
    # states specific to SNOWTRAN3D
    read!(joinpath(folder, "states_in_hiwl.bin"), histowet)
    read!(joinpath(folder, "states_in_sblt.bin"), dSWE_tot_subl)
    read!(joinpath(folder, "states_in_sltt.bin"), dSWE_tot_salt)
    read!(joinpath(folder, "states_in_sspt.bin"), dSWE_tot_susp)
  end

  if (fsm.SNSLID == 1)
    # states specific to SnowSlide
    read!(joinpath(folder, "landuse_slope.bin"), slope)
    read!(joinpath(folder, "landuse_shd.bin"), Shd)
    read!(joinpath(folder, "states_in_sldt.bin"), dSWE_tot_slide)
    read!(joinpath(folder, "states_in_idem.bin"), index_sorted_dem)
  end

  # Tuned snow surface properties

  if fsm.OSHDTN == 0

    fsm.adm = 100
    fsm.adc[:,:] = 1000
    fsm.afs[:,:] = fsm.asmx
    if (fsm.TILE == "glacier" || ((fsm.SNTRAN == 1 || fsm.SNSLID == 1) && fsm.glacierfrac(i,j) > eps(Tf)))
      fsm.z0_snow[:,:] = 0.0009
    else
      fsm.z0_snow[:,:] = fsm.z0sn
    end

  else

    fsm.adm = 130

    for j = 1:fsm.Ny
      for i = 1:fsm.Nx
      
        # Elevation-dependent tuning of cold snow albedo decay time
        if (fsm.dem[i,j] >= 2300)
          fsm.adc[i,j]  = 6000
        elseif (fsm.dem[i,j] <= 1500)
          fsm.adc[i,j] = 3000
        else
          fsm.adc[i,j] = 6000 + (2300 - fsm.dem[i,j]) / (2300 - 1500) * (3000 - 6000)
        end
        
        # Fresh snow albedo is now constant (previously elevation-dependent)
        fsm.afs[i,j] = fsm.asmx
        
        # Elevation-dependent tuning of snow roughness length
        if (fsm.TILE == "glacier" || ((fsm.SNTRAN == 1 || fsm.SNSLID == 1) && fsm.glacierfrac[i,j] > eps(Tf)))
          fsm.z0_snow[i,j] = 0.0009
        elseif (fsm.TILE == "forest")
          fsm.z0_snow[i,j] = fsm.z0sn
        else
          if (fsm.dem[i,j] >= 2300)
            fsm.z0_snow[i,j] = 0.01
          elseif (fsm.dem[i,j] >= 1500)
            fsm.z0_snow[i,j] = 0.2 + (fsm.dem[i,j] - 1500) / (2300 - 1500) * (0.01 - 0.2)
          else
            fsm.z0_snow[i,j] = 0.2
          end
        end
        
      end
    end

  end

  return fsm

end
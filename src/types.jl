@with_kw mutable struct FSM{Tfloat, Tint}

  # Value of undefined variables

  undef = 1.0e6
  
  # Maximum snow and soil layer thicknesses (m)

  Dzsnow::Vector{Tfloat} = [0.1, 0.2, 0.4]
  Dzsoil::Vector{Tfloat} = [0.1, 0.2, 0.4, 0.8]

  # Number of snow and soil layers

  Nsmax::Tint = length(Dzsnow)
  Nsoil::Tint = length(Dzsoil)
  
  # Domain size

  Nx::Tint = 1
  Ny::Tint = 1

  # Driving data

  dt::Tfloat = 3600                                                   # Time step (s)
  zT::Tfloat = 10                                                     # Temperature measurement height (m)
  zU::Tfloat = 10                                                     # Wind speed measurement height (m)
  zRH::Tfloat = 10                                                    # Relative humidity measurement height (m)

  # Model configuration

  ALBEDO::Tint = 2                                                    # Snow albedo (0, 1, 2)
  CANMOD::Tint = 0                                                    # Forest canopy (0, 1)
  CONDCT::Tint = 1                                                    # Snow thermal conductivity (0, 1)
  DENSTY::Tint = 3                                                    # Snow density (0, 1, 2, 3)
  EXCHNG::Tint = 1                                                    # Turbulent exchange (0, 1)
  HYDROL::Tint = 2                                                    # Snow hydraulics (0, 1, 2)
  SNFRAC::Tint = 3                                                    # Snow cover fraction (0, 1, 2, 3, 4)
  RADSBG::Tint = 0                                                    # Subgrid radiation param (0, 1)
  ZOFFST::Tint = 0                                                    # Measurement height offset (0, 1)
  OSHDTN::Tint = 1                                                    # OSHD-specific tuning options (0, 1)
  ALRADT::Tint = 0                                                    # Activate tuning of albedo decay as function of incoming direct SWR (0, 1)
  SNTRAN::Tint = 0                                                    # Snow transport (0, 1)
  SNSLID::Tint = 0                                                    # Snow slides (0, 1)
  SNOLAY::Tint = 0                                                    # Density-dependent layering (0, 1)
  CHECKS::Tint = 0                                                    # Check state variables at every time step (0, 1, 2)
  HN_ON::Bool = false        # TODO remove?                           # Activate the new snow model
  FOR_HN::Bool = true        # TODO remove?                           # Write 18h states for the HN model
  Z0PERT::Bool = false       # TODO remove?                           # Activate z0 perturbations
  WCPERT::Bool = false       # TODO remove?                           # Activate liquid water capacity perturbations
  FSPERT::Bool = false       # TODO remove?                           # Activate fresh snow density perturbations
  ALPERT::Bool = false       # TODO remove?                           # Activate albedo perturbations
  SLPERT::Bool = false       # TODO remove?                           # Activate settling perturbations

  # Tile options

  TILE::String = "open"                                               # Tile type
  tthresh::Tfloat = 0.1                                               # Tile threshold

  # Numerical solution parameters

  Nitr = 4                                                            # Number of iterations for surface energy balance

  # Defaults for canopy parameters

  avg0::Tfloat = 0.1                                                  # Snow-free vegetation albedo
  avgs::Tfloat = 0.4                                                  # Snow-covered vegetation albedo
  cden::Tfloat = 0.004                                                # Dense canopy turbulent transfer coefficient
  cvai::Tfloat = 4.4                                                  # Canopy snow capacity per unit VAI (kg/m^2)
  cveg::Tfloat = 20                                                   # Vegetation turbulent transfer coefficient
  Gcn1::Tfloat = 0.5                                                  # Leaf angle distribution parameter
  Gcn2::Tfloat = 0                                                    # Leaf angle distribution parameter
  gsnf::Tfloat = 0                                                    # Snow-free vegetation moisture conductance (m/s)
  kdif::Tfloat = 0.5                                                  # Diffuse radiation extinction coefficient
  kveg::Tfloat = 1                                                    # Canopy cover coefficient
  rchd::Tfloat = 0.67                                                 # Ratio of displacement height to canopy height
  rchz::Tfloat = 0.2                                                  # Ratio of roughness length to canopy height
  tcnc::Tfloat = 3600*240                                             # Canopy unloading time scale for cold snow (s)
  tcnm::Tfloat = 3600*48                                              # Canopy unloading time scale for melting snow (s)

  # Defaults for snow parameters

  a_eta::Tfloat = 0.1                                                 # Temperature factor for Crocus B92 compaction (K^-1)
  asmx::Tfloat = 0.86                                                 # Maximum albedo for fresh snow
  asmn::Tfloat = 0.6                                                  # Minimum albedo for melting snow
  b_eta::Tfloat = 0.023                                               # First density factor for Crocus B92 compaction (m^3/kg)
  bthr::Tfloat = 2                                                    # Snow thermal conductivity exponent
  c_eta::Tfloat = 250                                                 # Second density factor for Crocus B92 compaction (kg/m^3)
  eta0::Tfloat = 3.7e7                                                # Reference snow viscosity (Pa s)
  eta1::Tfloat = 7.62237e6                                            # Reference snow viscosity for Crocus B92 compaction (Pa s)
  hfsn::Tfloat = 0.1                                                  # Snowcover fraction depth scale (m)
  kfix::Tfloat = 0.24                                                 # Fixed thermal conductivity of snow (W/m/K)
  rho0::Tfloat = 300                                                  # Fixed snow density (kg/m^3)
  rhob::Tfloat = 6                                                    # Temperature factor in fresh snow density (kg/m^3/K)
  rhoc::Tfloat = 26                                                   # Wind factor in fresh snow density (kg s^0.5/m^3.5)
  rhof::Tfloat = 109                                                  # Fresh snow density (kg/m^3)
  rhos_min::Tfloat = 50                                               # Minimum snow density (kg/m^3)
  rhos_max::Tfloat = 750                                              # Maximum snow density (kg/m^3)
  rcld::Tfloat = 300                                                  # Maximum density for cold snow (kg/m^3)
  rgr0::Tfloat = 5e-5                                                 # Fresh snow grain radius (m)
  rmlt::Tfloat = 500                                                  # Maximum density for melting snow (kg/m^3)
  Salb::Tfloat = 10                                                   # Albedo decay constant (kg/m^2)
  snda::Tfloat = 2.8e-6                                               # Thermal metamorphism parameter (1/s)
  Talb::Tfloat = -2                                                   # Albedo decay temperature threshold (C)
  tcld::Tfloat = 3600*1000                                            # Cold snow albedo decay time scale (s)
  tmlt::Tfloat = 3600*100                                             # Melting snow albedo decay time scale (s)
  trho::Tfloat = 3600*200                                             # Snow compaction time scale (s)
  Wirr::Tfloat = 0.03                                                 # Irreducible liquid water content of snow
  z0sn::Tfloat = 0.002                                                # Snow roughness length (m)
  Sfmin::Tfloat = 10                                                  # Minimum 24h snowfall to refresh albedo (kg/m^2)

  # Defaults for ground surface parameters

  bstb::Tfloat = 5                                                    # Atmospheric stability parameter
  gsat::Tfloat = 0.01                                                 # Surface conductance for saturated soil (m/s)

  # Defaults for additional forest snow process parametrization

  adfs::Tfloat = 3                                                    # Snow albedo adjustment dependent on SWR
  adfl::Tfloat = 2                                                    # Snow albedo adjustment dependent on LWR
  fsar::Tfloat = 0.1                                                  # Snow albedo adjustment range dependent on vegetation fraction
  psf::Tfloat  = 1                                                    # Scaling factor for solid precipitation (within forest stand, at min CC)
  psr::Tfloat  = 0.1                                                  # Range of solid precipitation (within forest stand, spread min-max CC)
  wcan::Tfloat = 2.5                                                  # Parameter of exponential wind profile
  zsub::Tfloat = 2                                                    # Sub-canopy reference height (m)
  zgf::Tfloat = 1                                                     # Roughness length adjustment factor depending on vegetation fraction
  zgr::Tfloat = 0                                                     # Roughness length adjustment range depending on vegetation fraction
  khcf::Tfloat = 3                                                    # Diffusivity adjustment for canopy effects (Finnigan 2000)

  # Surface parameters

  adm::Tfloat = -999  # TODO defaults?                                # Melting snow albedo decay time (h)
  adc::Array{Tfloat,2} = -999*ones(Nx,Ny)  # TODO defaults?           # Cold snow albedo decay time (h)
  afs::Array{Tfloat,2} = -999*ones(Nx,Ny)  # TODO defaults?           # Maximum albedo for fresh snow
  z0_snow::Array{Tfloat,2} = -999*ones(Nx,Ny)  # TODO defaults?       # Roughness length of snow (m)

  # Surface properties

  alb0::Array{Tfloat,2} = 0.2*ones(Nx,Ny)                             # Snow-free ground albedo
  z0sf::Array{Tfloat, 2} = 0.2*ones(Nx,Ny)                            # Snow-free roughness length (m)
  fcly::Array{Tfloat, 2} = 0.3*ones(Nx,Ny)                            # Soil clay fraction
  fsnd::Array{Tfloat, 2} = 0.6*ones(Nx,Ny)                            # Soil sand fraction

  # Canopy parameters

  canh::Array{Tfloat, 2} = Tfloat(undef)*ones(Nx,Ny)                  # Canopy heat capacity (J/K/m^2)
  fsky::Array{Tfloat, 2} = Tfloat(undef)*ones(Nx,Ny)                  # Sky view fraction
  fveg::Array{Tfloat, 2} = Tfloat(undef)*ones(Nx,Ny)                  # Canopy cover fraction
  fves::Array{Tfloat, 2} = Tfloat(undef)*ones(Nx,Ny)                  # Stand-scale canopy cover fraction
  hcan::Array{Tfloat, 2} = Tfloat(undef)*ones(Nx,Ny)                  # Canopy height (m)
  lai::Array{Tfloat, 2} = Tfloat(undef)*ones(Nx,Ny)                   # Leaf area index
  pmultf::Array{Tfloat, 2} = Tfloat(undef)*ones(Nx,Ny)                # Precip multiplier to revert precip correction applied to open area
  scap::Array{Tfloat, 2} = Tfloat(undef)*ones(Nx,Ny)                  # Canopy snow capacity (kg/m^2)
  trcn::Array{Tfloat, 2} = Tfloat(undef)*ones(Nx,Ny)                  # Canopy transmissivity
  VAI::Array{Tfloat, 2} = Tfloat(undef)*ones(Nx,Ny)                   # Vegetation area index
  vfhp::Array{Tfloat, 2} = Tfloat(undef)*ones(Nx,Ny)                  # Hemispherical sky-view fraction including canopy

  # Terrain properties

  slopemu::Array{Tfloat, 2} = Tfloat(undef)*ones(Nx,Ny)               # Slope parameter
  xi::Array{Tfloat, 2} = Tfloat(undef)*ones(Nx,Ny)                    # Terrain correlation length
  Ld::Array{Tfloat, 2} = Tfloat(undef)*ones(Nx,Ny)                    # Grid cell size or domain size (m)
  lat::Array{Tfloat, 2} = Tfloat(undef)*ones(Nx,Ny)                   # Latitude of each grid cell (center?)
  lon::Array{Tfloat, 2} = Tfloat(undef)*ones(Nx,Ny)                   # Longitude of each grid cell (center?)
  dem::Array{Tfloat, 2} = Tfloat(undef)*ones(Nx,Ny)                   # Grid elevation (m)
  tilefrac::Array{Tfloat,2} = Tfloat(undef)*ones(Nx,Ny)               # Tile fraction
  glacierfrac::Array{Tfloat,2} = Tfloat(undef)*ones(Nx,Ny)            # Glacier flag

  # Derived soil parameters

  b::Array{Tfloat, 2} = zeros(Nx,Ny)                                  # Clapp-Hornberger exponent
  hcap_soil::Array{Tfloat, 2} = zeros(Nx,Ny)                          # Volumetric heat capacity of dry soil (J/K/m^3)
  hcon_soil::Array{Tfloat, 2} = zeros(Nx,Ny)                          # Thermal conductivity of dry soil (W/m/K)
  sathh::Array{Tfloat, 2} = zeros(Nx,Ny)                              # Saturated soil water pressure (m)
  Vsat::Array{Tfloat, 2} = zeros(Nx,Ny)                               # Volumetric soil moisture at saturation
  Vcrit::Array{Tfloat, 2} = zeros(Nx,Ny)                              # Volumetric soil moisture at critical point

  # State variables

  albs::Array{Tfloat, 2} = Tfloat(0.85)*ones(Nx,Ny)                   # Snow albedo
  Ds::Array{Tfloat,3} = zeros(Nsmax,Nx,Ny)                            # Snow layer thicknesses (m)
  Nsnow::Array{Tint,2} = zeros(Tint,Nx,Ny)                            # Number of snow layers
  Qcan::Array{Tfloat, 2} = Tfloat(undef)*ones(Nx,Ny)                  # Canopy air space humidity
  rgrn::Array{Tfloat,3} = Tfloat(undef)*ones(Nsmax,Nx,Ny)             # Snow layer grain radius (m)
  histowet::Array{Tfloat,3} = Tfloat(undef)*ones(Nsmax,Nx,Ny)         # Historical variable for past wetting of a layer (0-1)
  Sice::Array{Tfloat,3} = zeros(Nsmax,Nx,Ny)                          # Ice content of snow layers (kg/m^2)
  Sliq::Array{Tfloat,3} = zeros(Nsmax,Nx,Ny)                          # Liquid content of snow layers (kg/m^2)
  Sveg::Array{Tfloat, 2} = Tfloat(undef)*ones(Nx,Ny)                  # Snow mass on vegetation (kg/m^2)
  Tcan::Array{Tfloat, 2} = Tfloat(undef)*ones(Nx,Ny)                  # Canopy air space temperature (K)
  theta::Array{Tfloat,3} = zeros(Nsoil,Nx,Ny)                         # Volumetric moisture content of soil layers
  Tsnow::Array{Tfloat,3} = Tfloat(273.15)*ones(Nsmax,Nx,Ny)           # Snow layer temperatures (K)
  Tsoil::Array{Tfloat,3} = Tfloat(285)*ones(Nsoil,Nx,Ny)              # Soil layer temperatures (K)
  Tsrf::Array{Tfloat, 2} = Tfloat(285)*ones(Nx,Ny)                    # Surface skin temperature (K)
  fsnow::Array{Tfloat, 2} = zeros(Nx,Ny)                              # Snow cover fraction terrain
  Tveg::Array{Tfloat, 2} = Tfloat(undef)*ones(Nx,Ny)                  # Vegetation temperature (K)
  snowdepthmin::Array{Tfloat, 2} = zeros(Nx,Ny)                       # Minimum snow depth at time step of swemin (m)
  snowdepthmax::Array{Tfloat, 2} = zeros(Nx,Ny)                       # Maximum snow depth at time stemp of swemax (m)
  snowdepthhist::Array{Tfloat,3} = zeros(14,Nx,Ny)                    # History of snow depth during last 14 days (m). Most recent entries first
  swemin::Array{Tfloat, 2} = zeros(Nx,Ny)                             # Minimum SWE during the season (m)
  swemax::Array{Tfloat, 2} = zeros(Nx,Ny)                             # Maximum SWE during the season (m)
  swehist::Array{Tfloat,3} = zeros(14,Nx,Ny)                          # History of SWE during last 14 days (kg/m^2). Most recent entries first
  fsky_terr::Array{Tfloat, 2} = Tfloat(undef)*ones(Nx,Ny)             # Sky view fraction terrain

  # Radiation - temporary arrays

  alb::Array{Tfloat,2} = zeros(Nx,Ny)
  asrf_out::Array{Tfloat,2} = zeros(Nx,Ny)
  Sdirt::Array{Tfloat,2} = zeros(Nx,Ny)
  Sdift::Array{Tfloat,2} = zeros(Nx,Ny)
  SWveg::Array{Tfloat,2} = zeros(Nx,Ny)
  SWsrf::Array{Tfloat,2} = zeros(Nx,Ny)
  SWsci::Array{Tfloat,2} = zeros(Nx,Ny)
  LWt::Array{Tfloat,2} = zeros(Nx,Ny)
  SWtopo_out::Array{Tfloat,2} = zeros(Nx,Ny)

  # Thermal - temporary arrays

  ksnow::Array{Tfloat,3} = zeros(Nsmax, Nx, Ny)
  csoil::Array{Tfloat,3} = zeros(Nsoil, Nx, Ny)
  ksoil::Array{Tfloat,3} = zeros(Nsoil, Nx, Ny)
  gs1::Array{Tfloat,2} = zeros(Nx, Ny)
  Ds1::Array{Tfloat,2} = zeros(Nx, Ny)
  Ts1::Array{Tfloat,2} = zeros(Nx, Ny)
  ks1::Array{Tfloat,2} = zeros(Nx, Ny)
  Tveg0::Array{Tfloat,2} = zeros(Nx, Ny)

  # Sfexch - temporary arrays

  KH::Array{Tfloat,2} = zeros(Nx, Ny)
  KHa::Array{Tfloat,2} = zeros(Nx, Ny)
  KHg::Array{Tfloat,2} = zeros(Nx, Ny)
  KHv::Array{Tfloat,2} = zeros(Nx, Ny)
  KWg::Array{Tfloat,2} = zeros(Nx, Ny)
  KWv::Array{Tfloat,2} = zeros(Nx, Ny)
  Usc::Array{Tfloat,2} = zeros(Nx, Ny)

  # Ebalsrf - temporary arrays

  dTs::Array{Tfloat, 2} = zeros(Nx,Ny)
  Esrf::Array{Tfloat,2} = zeros(Nx,Ny)
  Eveg::Array{Tfloat,2} = zeros(Nx,Ny)
  G::Array{Tfloat,2} = zeros(Nx,Ny)
  H::Array{Tfloat,2} = zeros(Nx,Ny)
  Hsrf::Array{Tfloat,2} = zeros(Nx,Ny)
  LE::Array{Tfloat,2} = zeros(Nx,Ny)
  LEsrf::Array{Tfloat,2} = zeros(Nx,Ny)
  LWsci::Array{Tfloat,2} = zeros(Nx,Ny)
  LWveg::Array{Tfloat,2} = zeros(Nx,Ny)
  Melt::Array{Tfloat,2} = zeros(Nx,Ny)
  Rnet::Array{Tfloat,2} = zeros(Nx,Ny)
  Rsrf::Array{Tfloat,2} = zeros(Nx,Ny)
  Ssub::Tfloat = zero(Tfloat)

  # Canpopy - temporary arrays

  intcpt::Array{Tfloat,2} = zeros(Nx, Ny)
  Sbveg::Array{Tfloat,2} = zeros(Nx, Ny)
  unload::Array{Tfloat, 2} = zeros(Nx, Ny)

  # Snow - temporary arrays

  Gsoil::Array{Tfloat,2} = zeros(Nx, Ny)
  Roff::Array{Tfloat,2} = zeros(Nx, Ny)
  meltflux_out::Array{Tfloat,2} = zeros(Nx, Ny)
  Sbsrf::Array{Tfloat,2} = zeros(Nx, Ny)
  Roff_bare::Array{Tfloat, 2} = zeros(Nx, Ny)
  Roff_snow::Array{Tfloat, 2} = zeros(Nx, Ny)
  fsnow_thres::Array{Tfloat, 2} = zeros(Nx, Ny)

  a::Vector{Tfloat} = zeros(Nsmax)
  bsnow::Vector{Tfloat} = zeros(Nsmax)
  c::Vector{Tfloat} = zeros(Nsmax)
  csnow::Vector{Tfloat} = zeros(Nsmax)
  dTssnow::Vector{Tfloat} = zeros(Nsmax)
  D::Vector{Tfloat} = zeros(Nsmax)
  E::Vector{Tfloat} = zeros(Nsmax)
  Gs::Vector{Tfloat} = zeros(Nsmax)
  rhs::Vector{Tfloat} = zeros(Nsmax)
  R::Vector{Tfloat} = zeros(Nsmax)
  S::Vector{Tfloat} = zeros(Nsmax)
  U::Vector{Tfloat} = zeros(Nsmax)
  W::Vector{Tfloat} = zeros(Nsmax)

  SWEbuffer::Vector{Tfloat} = zeros(15)
  snowdepthbuffer::Vector{Tfloat} = zeros(15)
  diffSWEbuffer::Vector{Tfloat} = zeros(14)

  # Soil - temporary vectors

  asoil::Vector{Tfloat} = zeros(Nsoil)
  bsoil::Vector{Tfloat} = zeros(Nsoil)
  cssoil::Vector{Tfloat} = zeros(Nsoil)
  dTssoil::Vector{Tfloat} = zeros(Nsoil)
  Gssoil::Vector{Tfloat} = zeros(Nsoil)
  rhssoil::Vector{Tfloat} = zeros(Nsoil)

  # Tridiag - temporary vectors

  gammasnow::Vector{Tfloat} = zeros(Nsmax)
  gammasoil::Vector{Tfloat} = zeros(Nsoil)

end

@with_kw mutable struct MET{Tfloat, Tint}
  
  # Domain size

  Nx::Tint = 1
  Ny::Tint = 1
  
  # Time variables

  year::Array{Tint, 2} = -9999*ones(1, 1)
  month::Array{Tint, 2} = -9999*ones(1, 1)
  day::Array{Tint, 2} = -9999*ones(1, 1)
  hour::Array{Tfloat, 2} = -9999*ones(1, 1)

  # Meteorological variables

  Sdir::Array{Tfloat, 2} = zeros(Nx, Ny)                              # Direct-beam shortwave radiation (W/m^2)
  Sdif::Array{Tfloat, 2} = zeros(Nx, Ny)                              # Diffuse shortwave radiation (W/m^2)
  Sdird::Array{Tfloat, 2} = zeros(Nx, Ny)                             # Direct-beam shortwave radiation, per horizontal surface area (W/m^2)
  LW::Array{Tfloat, 2} = zeros(Nx, Ny)                                # Incoming longwave radiation (W/m^2)
  Sf::Array{Tfloat, 2} = zeros(Nx, Ny)                                # Snowfall rate (kg/m^2/s)
  Rf::Array{Tfloat, 2} = zeros(Nx, Ny)                                # Rainfall rate (kg/m^2/s)
  Ta::Array{Tfloat, 2} = zeros(Nx, Ny)                                # Air temperature (K)
  RH::Array{Tfloat, 2} = zeros(Nx, Ny)                                # Relative humidity (%)
  Ua::Array{Tfloat, 2} = zeros(Nx, Ny)                                # Wind speed (m/s)
  Udir::Array{Tfloat, 2} = zeros(Nx, Ny)                              # Wind direction (degrees, clockwise from N)
  Ps::Array{Tfloat, 2} = zeros(Nx, Ny)                                # Surface pressure (Pa)
  Sf24h::Array{Tfloat, 2} = zeros(Nx, Ny)                             # Snowfall 24hr (kg/m^2)
  Tc::Array{Tfloat, 2} = zeros(Nx, Ny)                                # Canopy temperature (K)
  es::Array{Tfloat, 2} = zeros(Nx, Ny)                                # Saturation vapour pressure (Pa)
  Qa::Array{Tfloat, 2} = zeros(Nx, Ny)                                # Specific humidity (kg/kg)
  Tv::Array{Tfloat, 2} = ones(Nx, Ny)                                 # Time-varying transmissivity for direct SWR (-)

  # Snowfall tracking variable

  Sf_history::Array{Tfloat, 3} = zeros(Nx, Ny, 24)

end

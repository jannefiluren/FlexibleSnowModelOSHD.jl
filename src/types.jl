@with_kw mutable struct FSM{Tfloat, Tint}
  
  # Maximum snow and soil layer thicknesses
  Dzsnow::Vector{Tfloat} = [0.1, 0.2, 0.4]
  Dzsoil::Vector{Tfloat} = [0.1, 0.2, 0.4, 0.8]
  
  # Base variables
  Nx::Tint = 1
  Ny::Tint = 1
  Nsmax::Tint = length(Dzsnow)
  Nsoil::Tint = length(Dzsoil)

  # Driving data
  dt::Tfloat = 3600*1
  zT::Tfloat = 10
  zU::Tfloat = 10
  zRH::Tfloat = 10

  # Model configuration
  ALBEDO::Tint = 2
  CANMOD::Tint = 0
  CONDCT::Tint = 1
  DENSTY::Tint = 3
  EXCHNG::Tint = 1
  HYDROL::Tint = 2
  SNFRAC::Tint = 3
  RADSBG::Tint = 0
  ZOFFST::Tint = 0
  OSHDTN::Tint = 1
  ALRADT::Tint = 0
  SNTRAN::Tint = 0
  SNSLID::Tint = 0
  SNOLAY::Tint = 0
  CHECKS::Tint = 0       # TODO remove?
  HN_ON::Bool = false
  FOR_HN::Bool = true
  Z0PERT::Bool = false
  WCPERT::Bool = false
  FSPERT::Bool = false
  ALPERT::Bool = false
  SLPERT::Bool = false

  tthresh::Tfloat = 0.1    # TODO check if this is the best handling
  TILE::String = "open"    # TODO check if this is the best handling and type

  # Numerical solution parameters
  Nitr = 4

  # Defaults for canopy parameters
  avg0::Tfloat = 0.1
  avgs::Tfloat = 0.4
  cden::Tfloat = 0.004
  cvai::Tfloat = 4.4
  cveg::Tfloat = 20
  Gcn1::Tfloat = 0.5
  Gcn2::Tfloat = 0
  gsnf::Tfloat = 0
  kdif::Tfloat = 0.5
  kveg::Tfloat = 1
  rchd::Tfloat = 0.67
  rchz::Tfloat = 0.2          
  tcnc::Tfloat = 3600*240
  tcnm::Tfloat = 3600*48

  # Defaults for snow parameters
  a_eta::Tfloat = 0.1
  asmx::Tfloat = 0.86
  asmn::Tfloat = 0.6
  b_eta::Tfloat = 0.023
  bthr::Tfloat = 2
  c_eta::Tfloat = 250
  eta0::Tfloat = 3.7e7
  eta1::Tfloat = 7.62237e6
  hfsn::Tfloat = 0.1 
  kfix::Tfloat = 0.24
  rho0::Tfloat = 300
  rhob::Tfloat = 6
  rhoc::Tfloat = 26
  rhof::Tfloat = 109
  rhos_min::Tfloat = 50
  rhos_max::Tfloat = 750
  rcld::Tfloat = 300
  rgr0::Tfloat = 5e-5
  rmlt::Tfloat = 500
  Salb::Tfloat = 10
  snda::Tfloat = 2.8e-6
  Talb::Tfloat = -2
  tcld::Tfloat = 3600*1000
  tmlt::Tfloat = 3600*100
  trho::Tfloat = 3600*200
  Wirr::Tfloat = 0.03
  z0sn::Tfloat = 0.002
  Sfmin::Tfloat = 10

  # Defaults for ground surface parameters
  bstb::Tfloat = 5
  gsat::Tfloat = 0.01

  # Defaults for additional parameters required for forest snow process parametrization
  adfs::Tfloat = 3
  adfl::Tfloat = 2
  fsar::Tfloat = 0.1
  psf::Tfloat  = 1
  psr::Tfloat  = 0.1
  wcan::Tfloat = 2.5
  zsub::Tfloat = 2
  zgf::Tfloat = 1
  zgr::Tfloat = 0
  khcf::Tfloat = 3

  # Surface parameters
  adm::Tfloat = -999  # TODO: what to do with defaults here
  adc::Array{Tfloat,2} = -999*ones(Nx,Ny)  # TODO: what to do with defaults here
  afs::Array{Tfloat,2} = -999*ones(Nx,Ny)  # TODO: what to do with defaults here
  z0_snow::Array{Tfloat,2} = -999*ones(Nx,Ny)  # TODO: what to do with defaults here

  # Surface properties
  alb0::Array{Tfloat,2} = 0.2*ones(Nx,Ny)
  z0sf::Array{Tfloat, 2} = 0.2*ones(Nx,Ny)
  fcly::Array{Tfloat, 2} = 0.3*ones(Nx,Ny)
  fsnd::Array{Tfloat, 2} = 0.6*ones(Nx,Ny)

  # Canopy parameters
  canh::Array{Tfloat, 2} = zeros(Nx,Ny)   # todo may need to change
  fsky::Array{Tfloat, 2} = zeros(Nx,Ny)   # todo may need to change
  fveg::Array{Tfloat, 2} = zeros(Nx,Ny)   # todo may need to change
  fves::Array{Tfloat, 2} = zeros(Nx,Ny)   # todo may need to change
  hcan::Array{Tfloat, 2} = zeros(Nx,Ny)   # todo may need to change
  lai::Array{Tfloat, 2} = zeros(Nx,Ny)   # todo may need to change
  pmultf::Array{Tfloat, 2} = zeros(Nx,Ny)   # todo may need to change
  scap::Array{Tfloat, 2} = zeros(Nx,Ny)   # todo may need to change
  trcn::Array{Tfloat, 2} = zeros(Nx,Ny)   # todo may need to change
  VAI::Array{Tfloat, 2} = zeros(Nx,Ny)   # todo may need to change
  vfhp::Array{Tfloat, 2} = zeros(Nx,Ny)   # todo may need to change

  # Terrain properties
  slopemu::Array{Tfloat, 2} = zeros(Nx,Ny)
  xi::Array{Tfloat, 2} = zeros(Nx,Ny)
  Ld::Array{Tfloat, 2} = zeros(Nx,Ny)
  lat::Array{Tfloat, 2} = zeros(Nx,Ny)
  lon::Array{Tfloat, 2} = zeros(Nx,Ny)
  dem::Array{Tfloat, 2} = zeros(Nx,Ny)
  tilefrac::Array{Tfloat,2} = zeros(Nx,Ny)
  glacierfrac::Array{Tfloat,2} = zeros(Nx,Ny)

  # Derived soil parameters
  b::Array{Tfloat, 2} = zeros(Nx,Ny)
  hcap_soil::Array{Tfloat, 2} = zeros(Nx,Ny)
  hcon_soil::Array{Tfloat, 2} = zeros(Nx,Ny)
  sathh::Array{Tfloat, 2} = zeros(Nx,Ny)
  Vsat::Array{Tfloat, 2} = zeros(Nx,Ny)
  Vcrit::Array{Tfloat, 2} = zeros(Nx,Ny)

  # State variables
  albs::Array{Tfloat, 2} = 0.85.*ones(Nx,Ny)
  Ds::Array{Tfloat,3} = zeros(Nsmax,Nx,Ny)
  Nsnow::Array{Tint,2} = zeros(Nx,Ny)
  Qcan::Array{Tfloat, 2} = zeros(Nx,Ny)
  rgrn::Array{Tfloat,3} = zeros(Nsmax,Nx,Ny)
  histowet::Array{Tfloat,3} = zeros(Nsmax,Nx,Ny)
  Sice::Array{Tfloat,3} = zeros(Nsmax,Nx,Ny)
  Sliq::Array{Tfloat,3} = zeros(Nsmax,Nx,Ny)
  Sveg::Array{Tfloat, 2} = zeros(Nx,Ny)
  Tcan::Array{Tfloat, 2} = 273.15.*ones(Nx,Ny)
  theta::Array{Tfloat,3} = zeros(Nsoil,Nx,Ny)
  Tsnow::Array{Tfloat,3} = 273.15.*ones(Nsmax,Nx,Ny)
  Tsoil::Array{Tfloat,3} = 273.15.*ones(Nsoil,Nx,Ny)
  Tsrf::Array{Tfloat, 2} = zeros(Nx,Ny)
  fsnow::Array{Tfloat, 2} = zeros(Nx,Ny)
  Tveg::Array{Tfloat, 2} = 273.15.*ones(Nx,Ny)
  snowdepthmin::Array{Tfloat, 2} = zeros(Nx,Ny)
  snowdepthmax::Array{Tfloat, 2} = zeros(Nx,Ny)
  snowdepthhist::Array{Tfloat,3} = zeros(14,Nx,Ny)
  swemin::Array{Tfloat, 2} = zeros(Nx,Ny)
  swemax::Array{Tfloat, 2} = zeros(Nx,Ny)
  swehist::Array{Tfloat,3} = zeros(14,Nx,Ny)
  fsky_terr::Array{Tfloat, 2} = zeros(Nx,Ny)

  # Radiation - temporary arrays

  #afs::Array{Tfloat,2} = zeros(Nx, Ny)
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
  #base variables
  Nx::Tint = 1
  Ny::Tint = 1
  year::Array{Tint, 2} = -9999*ones(1, 1)
  month::Array{Tint, 2} = -9999*ones(1, 1)
  day::Array{Tint, 2} = -9999*ones(1, 1)
  hour::Array{Tfloat, 2} = -9999*ones(1, 1)

  Sdir::Array{Tfloat, 2} = zeros(Nx, Ny)
  Sdif::Array{Tfloat, 2} = zeros(Nx, Ny)
  Sdird::Array{Tfloat, 2} = zeros(Nx, Ny)
  LW::Array{Tfloat, 2} = zeros(Nx, Ny)
  Sf::Array{Tfloat, 2} = zeros(Nx, Ny)
  Rf::Array{Tfloat, 2} = zeros(Nx, Ny)
  Ta::Array{Tfloat, 2} = zeros(Nx, Ny)
  RH::Array{Tfloat, 2} = zeros(Nx, Ny)
  Ua::Array{Tfloat, 2} = zeros(Nx, Ny)
  Udir::Array{Tfloat, 2} = zeros(Nx, Ny)
  Ps::Array{Tfloat, 2} = zeros(Nx, Ny)
  Sf24h::Array{Tfloat, 2} = zeros(Nx, Ny)
  Tc::Array{Tfloat, 2} = zeros(Nx, Ny)
  es::Array{Tfloat, 2} = zeros(Nx, Ny)
  Qa::Array{Tfloat, 2} = zeros(Nx, Ny)
  Tv::Array{Tfloat, 2} = ones(Nx, Ny)

  #arrays to calculate Sf24h
  Sf_history::Array{Tfloat, 3} = zeros(Nx, Ny, 24)
end
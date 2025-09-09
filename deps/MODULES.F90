!-----------------------------------------------------------------------
! Physical constants
!-----------------------------------------------------------------------
module CONSTANTS
real, parameter :: &
  cp = 1005,         &! Specific heat capacity of air (J/K/kg)
  eps = 0.622,       &! Ratio of molecular weights of water and dry air
  e0 = 610.78,       &! Saturation vapour pressure at Tm (Pa)
  grav = 9.81,       &! Acceleration due to gravity (m/s^2)
  hcap_ice = 2100,   &! Specific heat capacity of ice (J/K/kg)
  hcap_wat = 4180,   &! Specific heat capacity of water (J/K/kg)
  hcon_air = 0.025,  &! Thermal conductivity of air (W/m/K)
  hcon_clay = 1.16,  &! Thermal conductivity of clay (W/m/K)
  hcon_ice = 2.24,   &! Thermal conducivity of ice (W/m/K)
  hcon_sand = 1.57,  &! Thermal conductivity of sand (W/m/K)
  hcon_wat = 0.56,   &! Thermal conductivity of water (W/m/K)
  I0 = 1367,         &! Solar constant (W/m^2)
  Lf = 0.334e6,      &! Latent heat of fusion (J/kg)
  Lv = 2.501e6,      &! Latent heat of vapourisation (J/kg)
  Ls = Lf + Lv,      &! Latent heat of sublimation (J/kg)
  pi = 3.14159,      &! pi
  Rair = 287,        &! Gas constant for air (J/K/kg)
  Rwat = 462,        &! Gas constant for water vapour (J/K/kg)
  Runi = 8313,       &! Universal gas constant (J/kmol/K)
  rho_ice = 917,     &! Density of ice (kg/m^3)
  rho_wat = 1000,    &! Density of water (kg/m^3)
  rho_air = 1.275,   &! Density of air (kg/m^3)
  sb = 5.67e-8,      &! Stefan-Boltzmann constant (W/m^2/K^4)
  em_snow = 0.99,    &! Emissivity snow for Stefan-Boltzmann
  em_soil = 0.90,    &! Emissivity soil for Stefan-Boltzmann
  Tm = 273.15,       &! Melting point (K)
  vkman = 0.4,       &! Von Karman constant
  xM = 18.01,        &! Molecular weight of water (kg/kmol)
  visc_air = 13.e-6, &! Kinematic viscosity of air (m^2/s)
  undef = 1.e+6       ! Initialization value for allocatables
integer, parameter :: &
  iundef = 1.e+6      ! Initialization value for integer allocatables
end module CONSTANTS

!-----------------------------------------------------------------------
! Parameters for SnowTran3D
!-----------------------------------------------------------------------
module PARAM_SNOWTRAN3D
real, parameter :: &

  ! Define whether the threshold surface shear velocity will be
  ! constant during the simulation (Utau_t_flag = 0.0), or whether
  ! it will evolve as a funtion of air temperature and wind speed
  ! (Utau_t_flag = 1.0).
  !!! (recommended default value: Utau_t_flag = 0.0)
  Utau_t_flag = 1.0,    &! Flag for variable threshold friction velocitiy (1.0) or constant (0.0)

  ! For the case of Utau_t_flag = 0.0, define what the threshold
  ! surface shear velocity (m/s) will be during the simulation.
  !!! (recommended default value: Utau_t_const = 0.25)
  Utau_t_const = 0.25,  &! Constant threshold friction velocity (m/s) used if Utau_t_flag = 0.0

  ! Define whether you want two-layer accounting for soft (movable)
  ! and hard (unmovable) layers (flag = 0.0 = off, flag = 1.0 = on).
  !!! (recommended default value: twolayer_flag = 1.0)
  twolayer_flag = 1.0,  &! Flag for soft/hard layers distinction

  ! Define whether you want the upwind boundaries to have zero
  ! incoming transport flux, or the equilibrium transport flux.
  ! bc_flag = 0.0 corresponds to a zero flux, and bc_flag = 1.0
  ! corresponds to the equilibrium flux.
  !!! (recommended default value: bc_flag = 0.0)
  bc_flag = 0.0,        &! Boundary condition flag

  ! Density of snow deposited by snowdrift
  !!! (recommended default value: rho_snow = 300.0)
  rho_snow = 300.0,     &! Constant snow density (kg/m^3)

  ! The blowby parameter is implemented to account for the erosion
  ! of the tops of deep snow accumulations.  It corrects a
  ! deficiency in the du*/dx* < 0 formulation.  It is a number that
  ! should range from 0 to 1.0, and represents the fraction of the
  ! upwind saltation flux that is transfered farther downwind into
  ! the next grid cell.  So, the bigger the number, the less
  ! peaked the drift accumulation profile is.  blowby = 0.0 is the
  ! original model.  The Tabler surfaces can be used to do the
  ! same kind of thing.
  !!! (recommended default value: blowby = 0.0)
  blowby = 0.0           ! Fraction of the saltation flux transferred downwind
  
end module PARAM_SNOWTRAN3D

!-----------------------------------------------------------------------
! Constants for SnowTran3D
!-----------------------------------------------------------------------
module CONSTANTS_SNOWTRAN3D
real, parameter :: &

  ! Constants related to surface shear stress and saltation transport.
  fetch = 500.0,        &! Equilibrium fetch distance (m)
  xmu = 3.0,            &! Scaling constant for non-equilibrium saltation transport
  C_z = 0.12,           &! Coefficient 0.12 in Liston and Sturm (1998) eq. 5 p. 500
  h_const = 1.6,        &! Coefficient 1.6 in Liston and Sturm (1998) eq. 14 p. 501
  wind_min = 4.0,       &! Minimum wind speed to compute snow transport (m/s)
  
  ! Constants related to suspended snow profile.
  Up_const = 2.8,       &! Constant coefficient for calculation of U_p
  dz_susp = 0.20,       &! dz in suspension layer (m)
  ztop_susp = 2.0,      &! Height of the top of suspension layer (m)
  fall_vel = 0.3,       &! Particle-settling velocity, s in L&S 1998 eq. 12 (m/s)
  Ur_const = 0.5         ! Constant coefficient beta in L&S 1998, eq. 13
  
end module CONSTANTS_SNOWTRAN3D

!-----------------------------------------------------------------------
! Physical parameters for snow
!-----------------------------------------------------------------------
module PARAMETERS
real, parameter :: &
  rhos_min = 50.0,   &! Minimum snow density (kg/m^3)  
  rhos_max = 750.0    ! Maximum snow density (kg/m^3)
end module PARAMETERS

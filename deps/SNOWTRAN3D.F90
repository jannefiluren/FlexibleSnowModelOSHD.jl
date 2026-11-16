!-----------------------------------------------------------------------
! Snow transport by wind - Liston's scheme adapted to FSM 2.0
! References: Liston and Sturm (1998), Liston et al. (2007)
! Author: Glen E. Liston
! Adapted to FSM 2.0 by Louis Quéno
! Reference: Quéno et al. (2024)
! ATTENTION: y is the W->E axis, while x is the S->N
! South of (i,j): (i-1,j)
! North of (i,j): (i+1,j)
! West of (i,j): (i,j-1)
! East of (i,j): (i,j+1)
!-----------------------------------------------------------------------
subroutine SNOWTRAN3D(snowdepth0,Sice0,dSWE_salt,dSWE_susp,dSWE_subl)

use CONSTANTS, only: &
  pi,                      &! pi
  rho_wat,                 &! Density of water (kg/m^3)
  Tm                        ! Melting point (K)

use DRIVING, only: &  
  Ua,                      &! Wind speed (m/s)
  Udir                      ! Wind direction (degrees, clockwise from N)

use GRID, only: &
  Nx,Ny                     ! Grid dimensions

use STATE_VARIABLES, only: &
  Ds,                      &! Snow layer thicknesses (m)
  Nsnow,                   &! Number of snow layers
  fsnow,                   &! Snow cover fraction 
  Sice,                    &! Ice content of snow layers (kg/m^2)
  Sliq,                    &! Liquid content of snow layers (kg/m^2)
  Tsnow                     ! Snow layer temperatures (K)

use CONSTANTS_SNOWTRAN3D, only: &
  wind_min                  ! Minimum wind speed to compute snow transport (m/s)

implicit none

real, intent(inout) :: &
  snowdepth0(Nx,Ny),       &! Snow depth of snowdrift accumulation, averaged over the grid cell (m)
  Sice0(Nx,Ny),            &! Ice content of snowdrift accumulation(kg/m^2)
  dSWE_salt(Nx,Ny),        &! SWE change due to saltation (kg/m^2)
  dSWE_susp(Nx,Ny),        &! SWE change due to suspension (kg/m^2)
  dSWE_subl(Nx,Ny)          ! SWE change due to snowdrift sublimation (kg/m^2)

integer :: &
  index_wdir_E(Nx,2*Ny+1), &! Wind index array E
  index_wdir_W(Nx,2*Ny+1), &! Wind index array W
  index_wdir_N(Ny,2*Nx+1), &! Wind index array N
  index_wdir_S(Ny,2*Nx+1)   ! Wind index array S

logical :: &
  flag_blowing_snow         ! Blowing snow flag

real :: &
  Qsalt(Nx,Ny),            &! Saltation flux (kg/m/s)
  Qsalt_WE(Nx,Ny),         &! W->E (y) component of saltation flux (kg/m/s)
  Qsalt_SN(Nx,Ny),         &! S->N (x) component of saltation flux (kg/m/s)
  Qsusp_WE(Nx,Ny),         &! W->E (y) component of suspension flux (kg/m/s)
  Qsusp_SN(Nx,Ny),         &! S->N (x) component of suspension flux (kg/m/s)
  Qsubl(Nx,Ny),            &! Sublimation flux (kg/m^2/s)
  Ds_soft(Nx,Ny),          &! Soft snow thickness (m)
  snowthickness(Nx,Ny),    &! Depth of the snowpack in the snow covered part, i.e. not scaled by fsnow (m)
  wind_WE(Nx,Ny),          &! W->E (y) component of wind speed (m/s)
  wind_SN(Nx,Ny),          &! S->N (x) component of wind speed (m/s)
  Utau(Nx,Ny),             &! Friction velocity (m/s)
  Utau_t(Nx,Ny),           &! Threshold friction velocity (m/s)
  h_star(Nx,Ny),           &! Height of the saltation layer (m)
  z0(Nx,Ny)                 ! Surface roughness length (m)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize Ds_soft
Ds_soft(:,:) = 0.0

! Initialize all fluxes at 0
Qsalt(:,:) = 0.0
Qsalt_WE(:,:) = 0.0
Qsalt_SN(:,:) = 0.0
dSWE_salt(:,:) = 0.0
Qsusp_WE(:,:) = 0.0
Qsusp_SN(:,:) = 0.0
dSWE_susp(:,:) = 0.0
Qsubl(:,:) = 0.0
dSWE_subl(:,:) = 0.0

! Initialization of wind_WE and wind_SN
wind_WE = Ua * cos(-Udir * pi/180.0 - pi/2.0)
wind_SN = Ua * sin(-Udir * pi/180.0 - pi/2.0)

! Initialization of snowthickness
snowthickness = sum(Ds,dim=1)

! Update the thicknesses of the soft snow layer.
call COMPUTE_SOFT_SNOW(Ds_soft)

! Update the threshold friction velocity.
call COMPUTE_THRESHOLD_FRICTION_VELOCITY(Utau_t)

! Set the blowing snow flag to false until it is clear that we will have blowing snow.
flag_blowing_snow = .FALSE.

! If the wind speed is lower that some threshold, then don't
! need to to any of the snow transport computations.
if (maxval(Ua) >= wind_min) then

  ! Get the wind direction indexing arrays for this particular wind event (time step).
  call COMPUTE_WIND_DIRECTION_ARRAYS(index_wdir_E,index_wdir_W,index_wdir_N,index_wdir_S,wind_WE,wind_SN)

  ! Solve for Utau and z0 if snow is saltating, else solve assuming
  ! z0 is known from snow depth and/or veg type, and solve for Utau.
  call COMPUTE_FRICTION_VELOCITY(Utau,z0,h_star,snowthickness,flag_blowing_snow,Utau_t,Ds_soft)

  ! Update Ds_soft by removing too dense layers compared to U_tau
  call UPDATE_SOFT_SNOW(Utau,Ds_soft)

end if

! Compute snow transport fluxes and associated mass and depth changes
call FLUXES_AND_MASS_BALANCE(snowthickness,Ds_soft, &
                             index_wdir_E,index_wdir_W,index_wdir_N,index_wdir_S, &
                             flag_blowing_snow, &
                             wind_WE,wind_SN, &
                             Qsalt,Qsalt_WE,Qsalt_SN,dSWE_salt, &
                             Qsusp_WE,Qsusp_SN,dSWE_susp, &
                             Qsubl,dSWE_subl, &
                             snowdepth0,Sice0,Utau,Utau_t,z0,h_star)

end subroutine SNOWTRAN3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine COMPUTE_SOFT_SNOW(Ds_soft)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use STATE_VARIABLES, only: &
  Ds,                &! Snow layer thicknesses (m)
  Nsnow,             &! Number of snow layers
  histowet,          &! Historical variable for past wetting of a layer (0-1)
  Sliq                ! Liquid content of snow layers (kg/m^2)

implicit none

real, intent(inout) :: &
  Ds_soft(Nx,Ny)      ! Soft snow thickness (m)

integer :: &
  i,j,               &! x and y indexes
  k                   ! Layer counter

! New determination of soft snow based on layer wetting history and temperature
! Now histowet is used by default
do j = 1, Ny
  do i = 1, Nx

    k = 1
    do while (k <= Nsnow(i,j) .AND. Sliq(k,i,j) < epsilon(Sliq) .AND. histowet(k,i,j) < 0.5)
      Ds_soft(i,j) = Ds_soft(i,j) + Ds(k,i,j)
      k = k + 1
    end do

  end do
end do

end subroutine COMPUTE_SOFT_SNOW

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine UPDATE_SOFT_SNOW(Utau,Ds_soft)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use STATE_VARIABLES, only: &
  Ds,                &! Snow layer thicknesses (m)
  Nsnow,             &! Number of snow layers
  fsnow,             &! Snow cover fraction 
  Sice,              &! Ice content of snow layers (kg/m^2)
  Sliq                ! Liquid content of snow layers (kg/m^2)

implicit none

real, intent(in) :: &
  Utau(Nx,Ny)         ! Friction velocity (m/s)

real, intent(inout) :: &
  Ds_soft(Nx,Ny)      ! Soft snow thickness (m)

integer :: &
  i,j,               &! x and y indexes
  k                   ! Layer counter

logical :: &
  hard_layer_reached  ! A hard layer has been reached

real :: &
  Ds_soft_old,       &! Soft snow thickness before update (m)
  rho_layer,         &! Snow layer density
  Utau_t_layer        ! Threshold friction velocity of each layer (m/s)

do j = 1, Ny
  do i = 1, Nx

    Ds_soft_old = Ds_soft(i,j)
    Ds_soft(i,j) = 0.0
    hard_layer_reached = .FALSE.
    k = 1

    do while ((k <= Nsnow(i,j)) .AND. (Ds_soft_old - Ds_soft(i,j) > epsilon(Ds_soft)) .AND. (.NOT. hard_layer_reached))

      ! Calculate snow layer density
      rho_layer = (Sice(k,i,j) + Sliq(k,i,j)) / Ds(k,i,j) / fsnow(i,j)

      ! Calculate the snow threshold friction velocity.
      if (rho_layer <= 300.0) then
        Utau_t_layer = 0.10 * exp(0.003 * rho_layer)
      else
        Utau_t_layer = 0.005 * exp(0.013 * rho_layer)
      end if

      if (Utau(i,j) - Utau_t_layer > epsilon(Utau)) then
        Ds_soft(i,j) = Ds_soft(i,j) + Ds(k,i,j)
      else
        hard_layer_reached = .TRUE.
      end if

      k = k + 1

    end do

  end do
end do

end subroutine UPDATE_SOFT_SNOW

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine FLUXES_AND_MASS_BALANCE(snowthickness,Ds_soft, &
                                   index_wdir_E,index_wdir_W,index_wdir_N,index_wdir_S, &
                                   flag_blowing_snow, &
                                   wind_WE,wind_SN, &
                                   Qsalt,Qsalt_WE,Qsalt_SN,dSWE_salt, &
                                   Qsusp_WE,Qsusp_SN,dSWE_susp, &
                                   Qsubl,dSWE_subl, &
                                   snowdepth0,Sice0,Utau,Utau_t,z0,h_star)

use CONSTANTS, only: &
  rho_wat                   ! Density of water (kg/m^3)

use DRIVING, only: &
  dt                        ! Timestep (s)

use GRID, only: &
  Nx,Ny                     ! Grid dimensions

use LANDUSE, only: &
  tilefrac,                &! Grid cell tile fraction - here = open fraction since SNOWTRAN3D only runs in open tile
  glacierfrac               ! Glacier fraction

use PARAMMAPS, only: &
  veg_shd                   ! Vegetation snow holding depth (m)

use STATE_VARIABLES, only: &
  Ds,                      &! Snow layer thicknesses (m)
  fsnow,                   &! Snow cover fraction 
  Sice,                    &! Ice content of snow layers (kg/m^2)
  Sliq,                    &! Liquid content of snow layers (kg/m^2)
  dSWE_tot_subl,           &! Cumulated SWE change due to snowdrift sublimation (kg/m^2)
  dSWE_tot_salt,           &! Cumulated SWE change due to saltation (kg/m^2)
  dSWE_tot_susp             ! Cumulated SWE change due to suspension (kg/m^2)

use PARAM_SNOWTRAN3D, only: &
  flag_relax_suspension     ! Flag to relax suspension fluxes along the wind direction

implicit none

real, intent(in) :: &
  h_star(Nx,Ny),           &! Height of the saltation layer (m)
  z0(Nx,Ny),               &! Surface roughness length (m)
  wind_WE(Nx,Ny),          &! WE component of wind speed (m/s)
  wind_SN(Nx,Ny)            ! SN component of wind speed (m/s)

logical, intent(in) :: &
  flag_blowing_snow         ! Blowing snow flag

real, intent(inout) :: &
  Qsalt(Nx,Ny),            &! Saltation flux (kg/m/s)
  Qsalt_WE(Nx,Ny),         &! W->E (y) component of saltation flux (kg/m/s)
  Qsalt_SN(Nx,Ny),         &! S->N (x) component of saltation flux (kg/m/s)
  dSWE_salt(Nx,Ny),        &! SWE change due to saltation (kg/m^2)
  Qsusp_WE(Nx,Ny),         &! W->E (y) component of suspension flux (kg/m/s)
  Qsusp_SN(Nx,Ny),         &! S->N (x) component of suspension flux (kg/m/s)
  dSWE_susp(Nx,Ny),        &! SWE change due to suspension (kg/m^2)
  Qsubl(Nx,Ny),            &! Sublimation flux (kg/m^2/s)
  dSWE_subl(Nx,Ny),        &! SWE change due to snowdrift sublimation (kg/m^2)
  snowthickness(Nx,Ny),    &! Depth of the snowpack in the snow covered part, i.e. not scaled by fsnow (m)
  Ds_soft(Nx,Ny),          &! Soft snow thickness (m)
  snowdepth0(Nx,Ny),       &! Snow depth of snowdrift accumulation, averaged over the grid cell (m)
  Sice0(Nx,Ny),            &! Ice content of transported snow (kg/m^2)
  Utau(Nx,Ny),             &! Friction velocity (m/s)
  Utau_t(Nx,Ny)             ! Threshold friction velocity (m/s)

integer, intent(in) :: &
  index_wdir_E(Nx,2*Ny+1), &! Wind index array E
  index_wdir_W(Nx,2*Ny+1), &! Wind index array W
  index_wdir_N(Ny,2*Nx+1), &! Wind index array N
  index_wdir_S(Ny,2*Nx+1)   ! Wind index array S

integer :: &
  i,j                       ! Point counters

real :: &
  dSWE_subl_loss,          &! SWE loss due to sublimation (kg/m^2)
  dHS_subl_loss,           &! Snow depth loss due to sublimation (m)
  Ds_hard,                 &! Hard snow thickness (m)
  swe_loc,                 &! Total snowpack SWE (kg/m^2)
  snowdmin                  ! Minimum snow depth to allow transport (m)

real :: &
  openfrac(Nx,Ny),        &! Grid cell open fraction (including glaciers)
  Qsusp_eq_WE(Nx,Ny),     &! W->E (y) component of local-equilibrium suspension flux (kg/m/s)
  Qsusp_eq_SN(Nx,Ny),     &! S->N (x) component of local-equilibrium suspension flux (kg/m/s)
  Qsubl_susp(Nx,Ny),      &! Sublimation flux from the suspension layer (kg/m^2/s)
  Qsubl_salt(Nx,Ny),      &! Sublimation flux from the saltation layer (kg/m^2/s)
  Qsusp_eq,               &! Local-equilibrium suspension flux (kg/m/s)
  Qsusp,                  &! Suspension flux (kg/m/s)
  subl_ratio               ! Scaling of the suspension-layer sublimation: |Qsusp|/|Qsusp_eq|

! For snow transport, glaciers are considered open.
openfrac = tilefrac + glacierfrac

! If the blowing snow flag indicates wind transported snow
! somewhere within the domain (flag_blowing_snow = .TRUE.), run the saltation
! and suspension models.
if (flag_blowing_snow) then

  ! ---------------------------------------------------------------------------
  ! 1. SALTATION
  ! ---------------------------------------------------------------------------

  ! Solve the saltation flux.
  call SALTATION(Utau,Utau_t, &
                 index_wdir_E,index_wdir_W,index_wdir_N,index_wdir_S,wind_WE,wind_SN, &
                 Qsalt,Qsalt_WE,Qsalt_SN)

  ! Update SWE and HS changes due to saltation.
  ! Qsalt may be updated here due to limited erodible snow
  ! so it must be done before calculating suspension fluxes.
  call UPDATE_MASS_CHANGES(Qsalt_WE,Qsalt_SN,dSWE_salt, &
                           index_wdir_E,index_wdir_W,index_wdir_N,index_wdir_S, &
                           snowthickness,Ds_soft,snowdepth0,Sice0)

  ! Update Qsalt based on corrected WE and SN components, 
  ! even though this updated Qsalt is not needed later (to be clean).
  Qsalt = sqrt(Qsalt_WE**2 + Qsalt_SN**2)

  ! ---------------------------------------------------------------------------
  ! 2. SUSPENSION
  ! ---------------------------------------------------------------------------

  ! Solve the suspension flux (with potentially updated Qsalt).
  ! Qsusp_WE/Qsusp_SN hold the local-equilibrium (Kind diagnostic) flux.
  call SUSPENSION(Utau,z0,h_star,Utau_t,Qsalt_WE,Qsalt_SN, &
                  Qsusp_WE,Qsusp_SN,Qsubl_susp,Qsubl_salt)

  ! Save the local-equilibrium suspension fluxes.
  ! Even if no relaxation, they are used to scale the suspension-layer sublimation flux 
  ! because they were used to calculate Qsusp_subl before it may be adjusted later.
  Qsusp_eq_WE = Qsusp_WE
  Qsusp_eq_SN = Qsusp_SN

  if (flag_relax_suspension) then
    ! Relax the fluxes along the wind direction with fetch_susp
    call RELAX_SUSPENSION_FLUX(Qsusp_WE,Qsusp_SN, &
                              Qsusp_eq_WE,Qsusp_eq_SN, &
                              index_wdir_E,index_wdir_W,index_wdir_N,index_wdir_S)
  end if

  ! Update SWE and HS changes due to SUSPENSION.
  call UPDATE_MASS_CHANGES(Qsusp_WE,Qsusp_SN,dSWE_susp, &
                           index_wdir_E,index_wdir_W,index_wdir_N,index_wdir_S, &
                           snowthickness,Ds_soft,snowdepth0,Sice0)

  ! Adjust flux of sublimation in suspension because Qsusp may have been changed (either by relaxation or by limited erodible snow).
  ! The saltation-layer contribution is left unchanged because it was calculated after last changes of Qsalt.
  do j = 1, Ny
    do i = 1, Nx
      Qsusp_eq = sqrt(Qsusp_eq_WE(i,j)**2 + Qsusp_eq_SN(i,j)**2)
      Qsusp = sqrt(Qsusp_WE(i,j)**2 + Qsusp_SN(i,j)**2)
      if (Qsusp_eq > epsilon(Qsusp_eq)) then
        subl_ratio = Qsusp / Qsusp_eq
      else
        subl_ratio = 0.0
      end if
      Qsubl(i,j) = Qsubl_salt(i,j) + Qsubl_susp(i,j) * subl_ratio
    end do
  end do

  ! ---------------------------------------------------------------------------
  ! 3. SUBLIMATION
  ! ---------------------------------------------------------------------------

  do j = 1, Ny
    do i = 1, Nx

      ! Make adjustments for the case where there is no snow available
      ! on the ground (or captured within the vegetation) to be
      ! eroded.
      Ds_hard = snowthickness(i,j) - Ds_soft(i,j)
      snowdmin = max(veg_shd(i,j),Ds_hard)
      swe_loc = sum(Sice(:,i,j) + Sliq(:,i,j))

      ! Convert Qsubl to sublimated snow depth dHS_subl_loss.
      dSWE_subl(i,j) = Qsubl(i,j) * dt
      ! Qsubl represents the total tile sublimation since fsnow and openfrac scaling were already accounted for.
      ! To conserve mass, we need to divide by the open fraction of the pixel
      ! to compensate for the final tile combination multiplying by this fraction -> / openfrac(i,j)
      if (openfrac(i,j) >= 0.1) then
        dSWE_subl(i,j) = dSWE_subl(i,j) / openfrac(i,j)
      else
        dSWE_subl(i,j) = 0.0
      end if
      dSWE_subl_loss = - dSWE_subl(i,j)
      ! The following min should actually never be necessary, because handled by the snowdmin check later.
      ! Kept for clarity and safety.
      dSWE_subl_loss = min(dSWE_subl_loss,swe_loc)
      dSWE_subl(i,j) = - dSWE_subl_loss
      if (dSWE_subl_loss > epsilon(dSWE_subl_loss)) then
        call HS_FROM_SWE(dSWE_subl_loss,dHS_subl_loss,i,j)
      else
        dSWE_subl_loss = 0.0
        dHS_subl_loss = 0.0
        dSWE_subl(i,j) = 0.0
      end if

      if (snowthickness(i,j) > snowdmin .and. fsnow(i,j) > epsilon(fsnow)) then
        if (snowthickness(i,j) - dHS_subl_loss / fsnow(i,j) <= snowdmin) then
          dHS_subl_loss = (snowthickness(i,j) - snowdmin) * fsnow(i,j)
          if (dHS_subl_loss > epsilon(dHS_subl_loss)) then
            call SWE_FROM_HS(dHS_subl_loss,dSWE_subl_loss,i,j)
            dSWE_subl(i,j) = - dSWE_subl_loss
          else
            dHS_subl_loss = 0.0
            dSWE_subl_loss = 0.0
            dSWE_subl(i,j) = 0.0
          end if
        end if
      else
        dSWE_subl(i,j) = 0.0
        dSWE_subl_loss = 0.0
        dHS_subl_loss = 0.0
      end if

      ! Note that there is sublimation only on the open part of the pixel, but we don't need to weight it here
      ! since it is done in the final tile combination

      if (dSWE_subl_loss > epsilon(dSWE_subl_loss) .and. dHS_subl_loss > epsilon(dHS_subl_loss)) then
        call SNOW_ABLATION(dHS_subl_loss,dSWE_subl_loss,i,j)
        ! It is not necessary to update soft snow depth because it is not used anymore after this final step.
      end if

      ! Update the snow layer thicknesses
      snowthickness(i,j) = sum(Ds(:,i,j))

    end do
  end do

end if

! ---------------------------------------------------------------------------
! 4. Update cumulated SWE changes due to saltation, suspension and sublimation
! ---------------------------------------------------------------------------

dSWE_tot_salt= dSWE_tot_salt + dSWE_salt
dSWE_tot_susp = dSWE_tot_susp + dSWE_susp
dSWE_tot_subl = dSWE_tot_subl + dSWE_subl

end subroutine FLUXES_AND_MASS_BALANCE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SUSPENSION(Utau,z0,h_star,Utau_t,Qsalt_WE,Qsalt_SN, &
                      Qsusp_WE,Qsusp_SN,Qsubl_susp,Qsubl_salt)

use CONSTANTS, only: &
  vkman               ! Von Karman constant

use DRIVING, only: &
  Ta,                &! Air temperature (K)
  RH                  ! Relative humidity (%)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use CONSTANTS_SNOWTRAN3D, only: &
  Up_const,          &! Constant coefficient for calculation of U_p
  dz_susp,           &! dz in suspension layer
  ztop_susp,         &! Height of the top of suspension layer (m)
  fall_vel,          &! Particle-settling velocity, s in L&S 1998 eq. 12 (m/s)
  Ur_const            ! Constant coefficient beta in L&S 1998, eq. 13

implicit none

real, intent(in) :: &
  Qsalt_WE(Nx,Ny),      &! Saltation flux (kg/m/s)
  Qsalt_SN(Nx,Ny),      &! Saltation flux (kg/m/s)
  Utau_t(Nx,Ny),     &! Threshold friction velocity (m/s)
  Utau(Nx,Ny),       &! Friction velocity (m/s)
  h_star(Nx,Ny),     &! Height of the saltation layer (m)
  z0(Nx,Ny)           ! Surface roughness length (m)

real, intent(inout) :: &
  Qsusp_WE(Nx,Ny),    &! W->E (y) component of suspension flux (kg/m/s)
  Qsusp_SN(Nx,Ny),    &! S->N (x) component of suspension flux (kg/m/s)
  Qsubl_susp(Nx,Ny),  &! Sublimation flux from the suspension layer (kg/m^2/s)
  Qsubl_salt(Nx,Ny)    ! Sublimation flux from the saltation layer (kg/m^2/s)

integer :: &
  i,j,               &! x and y indexes
  nzsteps,           &! Number of vertical steps
  iz                  ! Vertical steps counter

real :: &
  conc,              &! Concentration of the suspended snow at height z (kg/m^3)
  conc_WE,           &! W->E (y) component of the concentration of the suspended snow at height z (kg/m^3)
  conc_SN,           &! S->N (x) component of the concentration of the suspended snow at height z (kg/m^3)
  conc_salt,         &! Saltation-layer reference-level mass concentration (kg/m^3)
  conc_salt_WE,      &! W->E (y) component of the saltation-layer reference-level mass concentration (kg/m^3)
  conc_salt_SN,      &! S->N (x) component of the saltation-layer reference-level mass concentration (kg/m^3)
  prd,               &! Auxiliary variable for calculations
  U_p,               &! Horizontal particle velocity within the saltation layer (m/s)
  U_r,               &! Wind speed at height h_star (m/s)
  phistar_Cr,        &! Left member of L&S 1998 eq. 13
  V_susp,            &! Sublimation loss rate coefficient in the suspension layer (s^-1)
  V_salt,            &! Sublimation loss rate coefficient in the saltation layer (s^-1)
  z                   ! Height (m)

! Compute the mass concentration of suspended snow according to
! Kind (1992).

do j = 1, Ny
  do i = 1, Nx

    conc_salt_WE = 0.0
    conc_salt_SN = 0.0
    conc_salt = 0.0

    if (Qsalt_WE(i,j) > epsilon(Qsalt_WE) .or. Qsalt_SN(i,j) > epsilon(Qsalt_SN)) then
      U_r = Utau(i,j)/vkman * log(h_star(i,j)/z0(i,j))
      phistar_Cr = Utau(i,j)/U_r * Ur_const
      prd = phistar_Cr * Utau(i,j) / fall_vel
      U_p = Up_const * Utau_t(i,j)

      ! Compute the concentration in the saltation layer (kg/m^3).
      ! Scaling by actual snow surface in open is already included in Qsalt
      ! hence inherited in conc_salt, then Qsusp and Qsubl
      conc_salt_WE = Qsalt_WE(i,j) / (h_star(i,j) * U_p)
      conc_salt_SN = Qsalt_SN(i,j) / (h_star(i,j) * U_p)
      conc_salt = sqrt(conc_salt_WE**2 + conc_salt_SN**2)

      nzsteps = int((ztop_susp - h_star(i,j)) / dz_susp)

      Qsusp_WE(i,j) = 0.0
      Qsusp_SN(i,j) = 0.0
      Qsubl_susp(i,j) = 0.0
      Qsubl_salt(i,j) = 0.0

      do iz = 1, nzsteps
        z = h_star(i,j) + 0.5 * dz_susp + real(iz - 1) * dz_susp

        ! Compute the concentration of the suspended snow at height z.
        conc_WE = conc_salt_WE * ((prd + 1.0) * &
               (z/h_star(i,j))**((-fall_vel)/(vkman*Utau(i,j))) - &
               prd)
        conc_WE = max(conc_WE,0.0)
        conc_SN = conc_salt_SN * ((prd + 1.0) * &
               (z/h_star(i,j))**((-fall_vel)/(vkman*Utau(i,j))) - &
               prd)
        conc_SN = max(conc_SN,0.0)
        conc = sqrt(conc_WE**2 + conc_SN**2)

        ! Only do the integration if the concentration is non-zero.
        if (conc_WE > epsilon(conc_WE)) then

          ! Perform the quadrature (summation), without the constants.
          if (z - z0(i,j) < epsilon(z)) then
            z = 1.2 * z0(i,j)
          end if
          Qsusp_WE(i,j) = Qsusp_WE(i,j) + conc_WE * log(z/z0(i,j)) * dz_susp

        end if

        ! Only do the integration if the concentration is non-zero.
        if (conc_SN > epsilon(conc_SN)) then

          ! Perform the quadrature (summation), without the constants.
          if (z - z0(i,j) < epsilon(z)) then
            z = 1.2 * z0(i,j)
          end if
          Qsusp_SN(i,j) = Qsusp_SN(i,j) + conc_SN * log(z/z0(i,j)) * dz_susp

        end if

        if (conc > epsilon(conc)) then
          ! Compute the sublimation due to suspension.
          call SUBLIMATION(z,RH(i,j),Ta(i,j),Utau(i,j), &
                            z0(i,j),V_susp,V_salt,Utau_t(i,j),.TRUE.)
          Qsubl_susp(i,j) = Qsubl_susp(i,j) + conc * V_susp * dz_susp
        end if

      end do

      ! Finish the quadratures.
      ! Include the constants for Qsusp.
      Qsusp_WE(i,j) = Utau(i,j) / vkman * Qsusp_WE(i,j)
      Qsusp_SN(i,j) = Utau(i,j) / vkman * Qsusp_SN(i,j)

      ! Include the sublimation contribution due to saltation.
      z = h_star(i,j) / 2.0
      call SUBLIMATION(z,RH(i,j),Ta(i,j),Utau(i,j), &
                       z0(i,j),V_susp,V_salt,Utau_t(i,j),.FALSE.)

      Qsubl_salt(i,j) = V_salt * conc_salt * h_star(i,j)

    else
      Qsusp_WE(i,j) = 0.0
      Qsusp_SN(i,j) = 0.0
      Qsubl_susp(i,j) = 0.0
      Qsubl_salt(i,j) = 0.0
    end if

  end do
end do

end subroutine SUSPENSION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine RELAX_SUSPENSION_FLUX(Qsusp_WE,Qsusp_SN,Qsusp_eq_WE,Qsusp_eq_SN, &
                                 index_wdir_E,index_wdir_W,index_wdir_N,index_wdir_S)

use GRID, only: &
  Nx,Ny                     ! Grid dimensions

use PARAM_SNOWTRAN3D, only: &
  flag_boundary_condition   ! Boundary condition flag

implicit none

real, intent(inout) :: &
  Qsusp_WE(Nx,Ny),         &! W->E (y) component of suspension flux (kg/m/s)
  Qsusp_SN(Nx,Ny)           ! S->N (x) component of suspension flux (kg/m/s)

real, intent(in) :: &
  Qsusp_eq_WE(Nx,Ny),      &! W->E (y) component of local-equilibrium suspension flux (kg/m/s)
  Qsusp_eq_SN(Nx,Ny)        ! S->N (x) component of local-equilibrium suspension flux (kg/m/s)

integer, intent(in) :: &
  index_wdir_E(Nx,2*Ny+1), &! Wind index array E
  index_wdir_W(Nx,2*Ny+1), &! Wind index array W
  index_wdir_N(Ny,2*Nx+1), &! Wind index array N
  index_wdir_S(Ny,2*Nx+1)   ! Wind index array S

! Addition to the original SnowTran-3D concept.
! Relaxation of the suspension flux toward its local equilibrium value, with advective memory along the wind direction.
! The diagnostic in SUSPENSION is a local equilibrium profile: it does not take into account any advective or non-steady effects.
! Without this relaxation, the suspended flux collapses within one grid cell behind crests and sheltered areas, concentrating the deposition.
! It mirrors the saltation treatment: the plume rebuilds toward the local equilibrium on a fetch_susp distance when the wind can carry more than it receives,
! and decays geometrically per grid cell when it receives more than the local wind can carry.

! ---------------------------------------------------------------------------
! 1. Upwind boundary condition
! ---------------------------------------------------------------------------

! Similar treatment as for the saltation fluxes: either zero incoming flux at the boundaries, 
! or steady-state (local equilibrium) incoming flux at the boundaries.
if (.NOT. flag_boundary_condition) then
  ! Zero incoming suspended flux at the boundaries.
  Qsusp_WE = 0.0
  Qsusp_SN = 0.0
else
  ! Steady-state (local equilibrium) incoming flux at the boundaries.
  Qsusp_WE = Qsusp_eq_WE
  Qsusp_SN = Qsusp_eq_SN
end if

! ---------------------------------------------------------------------------
! 2. Sweep each wind direction
!    wnd_dir: 1=West, 2=East, 3=South, 4=North
! ---------------------------------------------------------------------------

! Westerly winds (wnd_dir = 1)
call UPDATE_SUSPENSION_FLUX_PER_DIRECTION(Qsusp_WE,Qsusp_eq_WE,1,index_wdir_W,Nx,Ny)

! Easterly winds (wnd_dir = 2)
call UPDATE_SUSPENSION_FLUX_PER_DIRECTION(Qsusp_WE,Qsusp_eq_WE,2,index_wdir_E,Nx,Ny)

! Southerly winds (wnd_dir = 3)
call UPDATE_SUSPENSION_FLUX_PER_DIRECTION(Qsusp_SN,Qsusp_eq_SN,3,index_wdir_S,Ny,Nx)

! Northerly winds (wnd_dir = 4)
call UPDATE_SUSPENSION_FLUX_PER_DIRECTION(Qsusp_SN,Qsusp_eq_SN,4,index_wdir_N,Ny,Nx)

end subroutine RELAX_SUSPENSION_FLUX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine UPDATE_SUSPENSION_FLUX_PER_DIRECTION(Qsusp_dir,Qsusp_eq_dir,wnd_dir, &
                                               index_array,dim1_size,dim2_size)

use GRID, only: &
  Nx,Ny                                ! Grid dimensions

use LANDUSE, only: &
  cellsize,                           &! Grid cell size (m)
  tilefrac,                           &! Grid cell tile fraction - here = open fraction since SNOWSLIDE only runs in open tile
  glacierfrac                          ! Glacier fraction

use CONSTANTS_SNOWTRAN3D, only: &
  fetch_susp,                         &! Equilibrium fetch distance of the suspension layer (m)
  xmu_susp                             ! Scaling constant for non-equilibrium suspension transport

use FUNCTIONS_SNOWTRAN3D

implicit none

real, intent(inout) :: &
  Qsusp_dir(Nx,Ny)                     ! Suspension flux in the wind direction (kg/m/s)

real, intent(in) :: &
  Qsusp_eq_dir(Nx,Ny)                  ! Local-equilibrium suspension flux in the wind direction (kg/m/s)

integer, intent(in) :: &
  wnd_dir,                            &! Wind direction index (1=West, 2=East, 3=South, 4=North)
  dim1_size,                          &! Size of the first dimension
  dim2_size                            ! Size of the second dimension

integer, intent(in) :: &
  index_array(dim1_size,2*dim2_size+1) ! Wind index array

real :: &
  Qsusp_dir_in,                       &! Suspension flux in the wind direction at the upwind grid cell (kg/m/s)
  scale                                ! Scaling coefficient for the relaxation

integer :: &
  kstart, kend,                       &! Index boundaries
  incr,                               &! Loop increment (1 or -1 depending on wind direction)
  dim1, dim2,                         &! Indexes of first and second dimensions
  k,                                  &! Temporary index
  ic, jc,                             &! x and y indexes of the considered grid cell
  iuw, juw,                           &! x and y indexes of the upwind grid cell
  idw, jdw                             ! x and y indexes of the downwind grid cell

real :: &
  openfrac(Nx,Ny)                      ! Grid cell open fraction (including glaciers)

! Addition to the original SnowTran-3D concept.
! One directional sweep of the suspension-flux relaxation, mirroring UPDATE_SALTATION_FLUX_PER_DIRECTION.

! For snow transport, glaciers are considered open.
openfrac = tilefrac + glacierfrac

select case (wnd_dir)
case (1,3) ! Westerly and Southerly winds
  incr = 1
case (2,4) ! Easterly and Northerly winds
  incr = -1
end select

do dim1 = 1, dim1_size
  do k = 1, index_array(dim1,1)
    kstart = get_kstart(wnd_dir, dim1, k, index_array, dim1_size, dim2_size)
    kend   = get_kend(wnd_dir, dim1, k, index_array, dim1_size, dim2_size)

    do dim2 = kstart, kend, incr

      select case (wnd_dir)
      case (1,2) ! Westerly and Easterly winds
        ic = dim1
        jc = dim2
      case (3,4) ! Southerly and Northerly winds
        ic = dim2
        jc = dim1
      end select

      iuw = get_i_upwind(wnd_dir, ic)
      juw = get_j_upwind(wnd_dir, jc)
      idw = get_i_downwind(wnd_dir, ic)
      jdw = get_j_downwind(wnd_dir, jc)

      ! Scaling coefficient per grid cell (analogy of saltation applied to
      ! the suspension layer). Capped at 1.0 so the outgoing flux stays
      ! between the incoming flux and the local equilibrium flux.
      scale = min(1.0, xmu_susp * cellsize(ic,jc) / fetch_susp)

      ! Incoming flux from the upwind grid cell.
      Qsusp_dir_in = Qsusp_dir(iuw,juw)

      ! Symmetric relaxation toward the local equilibrium:
      ! - building zone (Q_in < Q_eq): rebuild on fetch_susp;
      ! - collapsing zone (Q_in > Q_eq): decay, implicit blowby = 1 - scale.
      Qsusp_dir(ic,jc) = Qsusp_dir_in + scale * (Qsusp_eq_dir(ic,jc) - Qsusp_dir_in)

      ! No flux from fully-forested cells (more than 90% forest)
      if (openfrac(ic,jc) < 0.1) then
        Qsusp_dir(ic,jc) = 0.0
      end if

      ! Downwind forest attenuates the flow
      if (idw >= 1 .and. idw <= Nx .and. &
          jdw >= 1 .and. jdw <= Ny) then
        ! If more than 90% forested downwind, the flux is fully blocked
        if (openfrac(idw,jdw) < 0.1) then
          Qsusp_dir(ic,jc) = 0.0
        else
          ! If partly forested downwind, attenuate the flux by the fraction of open area downwind
          Qsusp_dir(ic,jc) = Qsusp_dir(ic,jc) * openfrac(idw,jdw)
        end if
      end if

    end do
  end do
end do

end subroutine UPDATE_SUSPENSION_FLUX_PER_DIRECTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SALTATION(Utau,Utau_t, &
                     index_wdir_E,index_wdir_W,index_wdir_N,index_wdir_S,wind_WE,wind_SN, &
                     Qsalt,Qsalt_WE,Qsalt_SN)

use CONSTANTS, only: &
  grav,                    &! Acceleration due to gravity (m/s^2)
  rho_air                   ! Density of air (kg/m^3)

use GRID, only: &
  Nx,Ny                     ! Grid dimensions

use LANDUSE, only: &
  tilefrac,                &! Grid cell tile fraction - here = open fraction since SNOWSLIDE only runs in open tile
  glacierfrac               ! Glacier fraction

use STATE_VARIABLES, only: &
  fsnow                     ! Snow cover fraction 

use PARAM_SNOWTRAN3D, only: &
  flag_boundary_condition   ! Boundary condition flag

implicit none

real, intent(in) :: &
  wind_WE(Nx,Ny),          &! W->E (y) component of wind speed (m/s)
  wind_SN(Nx,Ny),          &! S->N (x) component of wind speed (m/s)
  Utau_t(Nx,Ny),           &! Threshold friction velocity (m/s)
  Utau(Nx,Ny)               ! Friction velocity (m/s)

real, intent(inout) :: &
  Qsalt(Nx,Ny),            &! Saltation flux (kg/m/s)
  Qsalt_WE(Nx,Ny),         &! W->E (y) component of saltation flux (kg/m/s)
  Qsalt_SN(Nx,Ny)           ! S->N (x) component of saltation flux (kg/m/s)

integer, intent(in) :: &
  index_wdir_E(Nx,2*Ny+1), &! Wind index array E
  index_wdir_W(Nx,2*Ny+1), &! Wind index array W
  index_wdir_N(Ny,2*Nx+1), &! Wind index array N
  index_wdir_S(Ny,2*Nx+1)   ! Wind index array S

integer :: &
  i,j                       ! x and y indexes

real :: &
  Qsalt_max(Nx,Ny),        &! Maximum possible saltation flux (kg/m/s)
  Qsalt_max_WE(Nx,Ny),     &! W->E (y) component of Qsalt_max (kg/m/s)
  Qsalt_max_SN(Nx,Ny),     &! S->N (x) component of Qsalt_max (kg/m/s)
  openfrac(Nx,Ny)           ! Grid cell open fraction (including glaciers)

! For snow transport, glaciers are considered open.
openfrac = tilefrac + glacierfrac

! ---------------------------------------------------------------------------
! 1. Maximum possible saltation flux (Eqn. 6, Liston & Sturm 1998)
! ---------------------------------------------------------------------------

! Compute the maximum possible saltation flux, assuming that
! an abundance of snow is available at the surface.
do j = 1, Ny
  do i = 1, Nx

    ! For a given wind speed, find Qsalt_max.
    Qsalt_max(i,j) = 0.68 * rho_air / grav * &
                     Utau_t(i,j) / Utau(i,j) * (Utau(i,j)**2 - Utau_t(i,j)**2)
    ! Saltation only occurs on the snow covered part (* fsnow) of the open fraction of the grid cell (* openfrac).
    Qsalt_max(i,j) = max(Qsalt_max(i,j)* fsnow(i,j) * openfrac(i,j),0.0)

    ! Now weight the max saltation flux for the u and v wind
    ! components, where the vector sum should equal Qsalt_max.
    Qsalt_max_WE(i,j) = Qsalt_max(i,j) * abs(wind_WE(i,j)) / &
                      sqrt(wind_WE(i,j)**2 + wind_SN(i,j)**2)
    Qsalt_max_SN(i,j) = Qsalt_max(i,j) * abs(wind_SN(i,j)) / &
                      sqrt(wind_WE(i,j)**2 + wind_SN(i,j)**2)
  end do
end do

! ---------------------------------------------------------------------------
! 2. Upwind boundary condition
! ---------------------------------------------------------------------------

! Define an upwind boundary condition.  If flag_boundary_condition = .TRUE. then it is
! assumed that the inflow saltation flux has reached steady state.
! If flag_boundary_condition = .FALSE. then the saltation flux is assumed to be zero.
! The boundary condition is implemented by initializing the arrays
! to Qsalt_max, and since upwind boundaries are not called in
! the Qsalt computation, they stay in effect for the future 
! accumulation/erosion computation.
if (.NOT. flag_boundary_condition) then
  ! Zero incoming flux at the boundaries.
  Qsalt_WE = 0.0
  Qsalt_SN = 0.0
else
  ! Steady-state (maximum) incoming flux at the boundaries.
  Qsalt_WE = Qsalt_max_WE
  Qsalt_SN = Qsalt_max_SN
end if

! ---------------------------------------------------------------------------
! 3. Sweep each wind direction
!    wnd_dir: 1=West, 2=East, 3=South, 4=North
! ---------------------------------------------------------------------------

! Westerly winds (wnd_dir = 1) — j ascending, upwind = j-1, downwind = j+1
call UPDATE_SALTATION_FLUX_PER_DIRECTION(Qsalt_WE, Qsalt_max_WE, Utau, 1, index_wdir_W, Nx, Ny)

! Easterly winds (wnd_dir = 2) — j descending, upwind = j+1, downwind = j-1
call UPDATE_SALTATION_FLUX_PER_DIRECTION(Qsalt_WE, Qsalt_max_WE, Utau, 2, index_wdir_E, Nx, Ny)

! Southerly winds (wnd_dir = 3) — i ascending, upwind = i-1, downwind = i+1
call UPDATE_SALTATION_FLUX_PER_DIRECTION(Qsalt_SN, Qsalt_max_SN, Utau, 3, index_wdir_S, Ny, Nx)

! Northerly winds (wnd_dir = 4) — i descending, upwind = i+1, downwind = i-1
call UPDATE_SALTATION_FLUX_PER_DIRECTION(Qsalt_SN, Qsalt_max_SN, Utau, 4, index_wdir_N, Ny, Nx)

! ---------------------------------------------------------------------------
! 4. Combine u and v components
! ---------------------------------------------------------------------------

! Combine the u and v components to yield the total saltation flux
! at each grid cell.
Qsalt = sqrt(Qsalt_WE**2 + Qsalt_SN**2)

end subroutine SALTATION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine UPDATE_SALTATION_FLUX_PER_DIRECTION(Qsalt_dir, Qsalt_max_dir, Utau, wnd_dir, index_array, dim1_size, dim2_size)

use GRID, only: &
  Nx,Ny                                ! Grid dimensions

use LANDUSE, only: &
  cellsize,                           &! Grid cell size (m)
  tilefrac,                           &! Grid cell tile fraction - here = open fraction since SNOWSLIDE only runs in open tile
  glacierfrac                          ! Glacier fraction

use PARAM_SNOWTRAN3D, only: &
  blowby                               ! Fraction of the saltation flux transferred downwind

use CONSTANTS_SNOWTRAN3D, only: &
  fetch,                              &! Equilibrium fetch distance (m)
  xmu                                  ! Scaling constant for non-equilibrium saltation transport

use FUNCTIONS_SNOWTRAN3D

implicit none

real, intent(inout) :: &
  Qsalt_dir(Nx,Ny)                     ! Saltation flux in the wind direction (kg/m/s)

real, intent(in) :: &
  Qsalt_max_dir(Nx,Ny),               &! Maximum possible saltation flux in the wind direction (kg/m/s)
  Utau(Nx,Ny)                          ! Friction velocity (m/s)

integer, intent(in) :: &
  wnd_dir,                            &! Wind direction index (1=West, 2=East, 3=South, 4=North)
  dim1_size,                          &! Size of the first dimension
  dim2_size                            ! Size of the second dimension

integer, intent(in) :: &
  index_array(dim1_size,2*dim2_size+1) ! Wind index array

real :: &
  dUtau,                              &! Utau difference (m/s) 
  Qsalt_dir_in,                       &! Saltation flux in the wind direction at upwind grid cell (kg/m/s)
  scale                                ! Scaling coefficient for Eqn. 9 in L&S 1998

integer :: &
  kstart, kend,                       &! Index boundaries
  incr,                               &! Loop increment (1 or -1 depending on wind direction)
  dim1, dim2,                         &! Indexes of first and second dimensions
  k,                                  &! Temporary index
  ic, jc,                             &! x and y indexes of the considered grid cell
  iuw, juw,                           &! x and y indexes of the upwind grid cell
  idw, jdw                             ! x and y indexes of the downwind grid cell

real :: &
  openfrac(Nx,Ny)                      ! Grid cell open fraction (including glaciers)

! For snow transport, glaciers are considered open.
openfrac = tilefrac + glacierfrac

select case (wnd_dir)
case (1,3) ! Westerly and Southerly winds
  incr = 1
case (2,4) ! Easterly and Northerly winds
  incr = -1
end select

do dim1 = 1, dim1_size
  do k = 1, index_array(dim1,1)
    kstart = get_kstart(wnd_dir, dim1, k, index_array, dim1_size, dim2_size)
    kend   = get_kend(wnd_dir, dim1, k, index_array, dim1_size, dim2_size)

    do dim2 = kstart, kend, incr

      select case (wnd_dir)
      case (1,2) ! Westerly and Easterly winds
        ic = dim1
        jc = dim2
      case (3,4) ! Southerly and Northerly winds
        ic = dim2
        jc = dim1
      end select

      iuw = get_i_upwind(wnd_dir, ic)
      juw = get_j_upwind(wnd_dir, jc)
      idw = get_i_downwind(wnd_dir, ic)
      jdw = get_j_downwind(wnd_dir, jc)

      ! Define the scaling coefficients for Eqn. 9 in L&S 1998. Don't
      ! let them be greater than 1.0 or you will make more snow than
      ! there was before.
      scale = min(1.0, xmu * cellsize(ic,jc) / fetch)

      ! Incoming flux from the upwind grid cell.
      Qsalt_dir_in = Qsalt_dir(iuw, juw)
      dUtau   = Utau(ic,jc) - Utau(iuw, juw)

      if (dUtau >= epsilon(dUtau)) then
        ! Wind increasing: relax toward equilibrium
        Qsalt_dir(ic,jc) = Qsalt_dir_in + scale * (Qsalt_max_dir(ic,jc) - Qsalt_dir_in)
      else
        ! Wind decreasing: carry flux forward or blowby
        if (Qsalt_dir_in < Qsalt_max_dir(ic,jc)) then
          Qsalt_dir(ic,jc) = Qsalt_dir_in
        else
          Qsalt_dir(ic,jc) = max(blowby * Qsalt_dir_in, Qsalt_max_dir(ic,jc))
        end if
      end if

      ! No flux from fully-forested cells (more than 90% forest)
      if (openfrac(ic,jc) < 0.1) then
        Qsalt_dir(ic,jc) = 0.0
      end if

      ! Downwind forest attenuates the flow
      if (idw >= 1 .and. idw <= Nx .and. &
          jdw >= 1 .and. jdw <= Ny) then
        ! If more than 90% forested downwind, the flux is fully blocked
        if (openfrac(idw,jdw) < 0.1) then
          Qsalt_dir(ic,jc) = 0.0
        else
          ! If partly forested downwind, attenuate the flux by the fraction of open area downwind
          Qsalt_dir(ic,jc) = Qsalt_dir(ic,jc) * openfrac(idw,jdw)
        end if
      end if

    end do
  end do
end do

end subroutine UPDATE_SALTATION_FLUX_PER_DIRECTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine COMPUTE_FRICTION_VELOCITY(Utau,z0,h_star,snowthickness,flag_blowing_snow,Utau_t,Ds_soft)

use CONSTANTS, only: &
  grav,                 &! Acceleration due to gravity (m/s^2)
  vkman                  ! Von Karman constant

use DRIVING, only: &
  Ua,                   &! Wind speed (m/s)
  zU                     ! Wind speed measurement height (m)

use GRID, only: &
  Nx,Ny                  ! Grid dimensions

use PARAMMAPS, only: &
  veg_shd,              &! Vegetation snow holding depth (m)
  z0sf,                 &! Snow-free roughness length (m)
  z0_snow                ! Roughness length of snow (m)

use STATE_VARIABLES, only: &
  fsnow                  ! Snow cover fraction 

use CONSTANTS_SNOWTRAN3D, only: &
  C_z,                  &! Coefficient 0.12 in Liston and Sturm (1998) eq. 5 p. 500
  h_const                ! Coefficient 1.6 in Liston and Sturm (1998) eq. 14 p. 501

implicit none

real, intent(in) :: &
  snowthickness(Nx,Ny), &! Depth of the snowpack in the snow covered part, i.e. not scaled by fsnow (m)
  Utau_t(Nx,Ny),        &! Threshold friction velocity (m/s)
  Ds_soft(Nx,Ny)         ! Soft snow thickness (m)

logical, intent(inout) :: &
  flag_blowing_snow      ! Blowing snow flag

real, intent(inout) :: &
  Utau(Nx,Ny)            ! Friction velocity (m/s)

real, intent(out) :: &  
  h_star(Nx,Ny),        &! Height of the saltation layer (m)
  z0(Nx,Ny)              ! Surface roughness length (m)

integer :: &
  i,j                    ! x and y indexes

logical :: &
  threshold_flag         ! Flag if u* is above threshold

real :: &
  threshold,            &! Threshold for u* (m/s)
  guess,                &! Initial guess for Utau
  Utautmp,              &! Temporary Utau variable (m/s)
  windtmp,              &! Temporary wind speed variable (m/s)
  wind_max,             &! Maximum wind speed (m/s)
  z0_tmp                 ! Temporary z0 variable (m)

! We are entering the subroutine with the blowing snow flag set to .FALSE.
! Then, if snow is found to blow in any
! domain grid cell, we set the flag to (flag_blowing_snow = .TRUE.).

! Build the Utau array.
guess = 0.1
do j = 1, Ny
  do i = 1, Nx

    ! Determine whether snow is saltating (this influences how Utau
    ! and z0 are computed).
    if (snowthickness(i,j) <= veg_shd(i,j) .or. Ds_soft(i,j) <= epsilon(Ds_soft)) then

      ! Saltation will not occur.
      z0(i,j) = fsnow(i,j) * z0_snow(i,j) + (1.0 - fsnow(i,j)) * z0sf(i,j) ! Eq. 4, LS 1998
      z0_tmp = min(0.25*zU,z0(i,j))
      Utau(i,j) = Ua(i,j) * vkman / log(zU/z0_tmp)
      h_star(i,j) = z0(i,j) * h_const / C_z

    else

      ! Saltation may occur.  Test for that possibility by assuming that
      ! saltation is present, solving for Utau and z0, and comparing
      ! whether Utau exceeds Utau_t.  If it does not, set z0 to that
      ! of snow and recompute Utau.

      ! To help insure that the iteration converges, set the minimum
      ! wind speed to be 1.0 m/s, and the maximum wind speed to be
      ! 30 m/s at 10-m height.
      windtmp = max(1.0,Ua(i,j))
      wind_max = 30.0 * log(zU/z0_snow(i,j))/log(10.0/z0_snow(i,j))
      windtmp = min(windtmp,wind_max) 

      ! For u* over 0.6, use the relation z0 = 0.00734 u* - 0.0022,
      ! instead of Equation (5) in Liston and Sturm (1998).  Note that
      ! for windspeeds greater than about 35 m/s this will have to be
      ! modified for the solution algorithm to converge (because the
      ! roughness length will start to be higher than the obs height!).
      threshold = 0.6/vkman * log(zU/0.0022)
      if (windtmp <= threshold) then
        threshold_flag = .FALSE.
      else
        threshold_flag = .TRUE.
      end if

      call FRICTION_VELOCITY_SOLVER(Utautmp,guess,windtmp,threshold_flag)

      if (Utautmp > Utau_t(i,j)) then

        ! We have saltation.
        Utau(i,j) = Utautmp
        z0(i,j) = C_z * Utau(i,j)**2 / (2.0 * grav)
        h_star(i,j) = h_const * Utau(i,j)**2 / (2.0 * grav)
        ! Set the blowing snow flag to .TRUE.
        flag_blowing_snow = .TRUE.

      else

        ! We do not have saltation, but the vegetation is covered by snow.
        ! Because we have determined that we do not have saltation, make
        ! sure Utau does not exceed Utau_t.
        z0(i,j) = z0_snow(i,j)
        Utau(i,j) = Ua(i,j) * vkman / log(zU/z0(i,j))
        Utau(i,j) = min(Utau(i,j),Utau_t(i,j))
        h_star(i,j) = z0(i,j) * h_const / C_z

      end if
    end if

  end do
end do

end subroutine COMPUTE_FRICTION_VELOCITY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine FRICTION_VELOCITY_SOLVER(xnew,guess,windtmp,threshold_flag)

use CONSTANTS, only: &
  grav,              &! Acceleration due to gravity (m/s^2)
  vkman               ! Von Karman constant

use DRIVING, only: &
  zU                  ! Wind speed measurement height (m)

use CONSTANTS_SNOWTRAN3D, only: &
  C_z                 ! Coefficient 0.12 in Liston and Sturm (1998) eq. 5 p. 500

implicit none

real, intent(in) :: &
  guess,             &! Initial guess for Utau (m/s)
  windtmp             ! Wind speed (m/s)

logical, intent(in) :: &
  threshold_flag      ! Flag if u* is above threshold

real, intent(inout) :: &  
  xnew                ! Solution for Utau (m/s)

integer :: &
  i,                 &! Point counter
  maxiter             ! Number of iterations

real :: &
  tol,               &! Convergence threshold
  old,               &! Temporary xnew
  fprime,            &! Temporary function 1
  funct               ! Temporary function 2

tol = 1.0e-3
maxiter = 20
old = guess

if (threshold_flag) then

  old = 0.6
  do i = 1, maxiter
    fprime = - 1.0 + windtmp * vkman * 0.00734 / (0.00734 * old - 0.0022) * &
            (log(zU) - log(0.00734 * old - 0.0022))**(-2)
    funct = - old + windtmp * vkman * &
           (log(zU) - log(0.00734 * old - 0.0022))**(-1)
    xnew = old - funct/fprime
    if (abs(xnew - old) < tol) return
    old = xnew
  end do

else

  do i = 1, maxiter
    fprime = - 1.0 + 2.0 / old * windtmp * vkman * &
            (log(zU) - log(C_z/(2.0*grav)) - 2.0*log(old))**(-2)
    funct = - old + windtmp * vkman * &
           (log(zU) - log(C_z/(2.0*grav)) - 2.0*log(old))**(-1)
    xnew = old - funct/fprime
    if (abs(xnew - old) < tol) return
    old = xnew
  end do

end if

write(*,*) 'max iteration exceeded when solving for Utau, Utau=',old

end subroutine FRICTION_VELOCITY_SOLVER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SUBLIMATION(z,RH,Ta,Utau,z0,V_susp,V_salt,Utau_t,flag_susp)

use CONSTANTS, only: &
  hcon_air,          &! Thermal conductivity of air (W/m/K)
  Ls,                &! Latent heat of sublimation (J/kg)
  pi,                &! pi
  Rair,              &! Gas constant for air (J/K/kg)
  Runi,              &! Universal gas constant (J/kmol/K)
  rho_ice,           &! Density of ice (kg/m^3)
  Tm,                &! Melting point (K)
  vkman,             &! Von Karman constant
  xM,                &! Molecular weight of water (kg/kmol)
  visc_air            ! Kinematic viscosity of air (m^2/s)

use DRIVING, only: &
  zRH                 ! Relative humidity measurement height (m)

implicit none

real, intent(in) :: &
  z,                 &! Height (m)
  RH,                &! Relative humidity (%)
  Ta,                &! Air temperature (K)
  Utau_t,            &! Threshold friction velocity (m/s)
  Utau,              &! Friction velocity (m/s)
  z0                  ! Surface roughness length (m)

logical, intent(in) :: &
  flag_susp           ! .TRUE.: sublimation due to suspension, .FALSE.: sublimation due to saltation

real, intent(out) :: &
  V_susp,            &! Sublimation loss rate coefficient in the suspension layer (s^-1)
  V_salt              ! Sublimation loss rate coefficient in the saltation layer (s^-1)

real :: &
  D,                 &! Diffusivity of water vapor in the atmosphere (m^2/s)
  rho_sat,           &! Saturation density of water vapour (kg/m^3)
  rh_offset,         &! Relative humidity offset (fraction 0->1)
  sigma,             &! Atmospheric undersaturation of water vapour with respect to ice at height z
  alpha,             &! Coefficient alpha defined in L&S (1998), eq. A-5 p. 515 
  rbar_r,            &! Mean radius of snow particles at height z (m)
  xmbar,             &! Mean particle mass at height z (kg)
  rbar,              &! Radius of a snow particle of mean particle mass xmbar at height z (m)
  u_z,               &! Wind speed at height z (m)
  x_r,               &! Fluctuating velocity component (m/s), defined in L&S (1998), eq. A-18 p. 516
  wbar,              &! Mean terminal fall velocity of suspended snow (m/s), L&S (1998), eq. A-17 p. 516
  V_r,               &! Ventilation velocity for turbulent suspension (m/s)
  V_rsalt,           &! Ventilation velocity for saltation (m/s)
  xN_r,              &! Particle Reynolds number Re
  xNu,               &! Nusselt number
  xSh,               &! Sherwood number
  tmp1,              &! Auxiliary variable 1
  tmp2,              &! Auxiliary variable 2
  top,               &! Top of eq. A-6, p. 516, L&S (1998)
  bottom              ! Bottom of eq. A-6, p. 516, L&S (1998)

D = 2.06e-5 * (Ta/273.)**(1.75)
! rho_sat = 0.622 * 10.0**(11.40 - 2353./Ta) / (Rair * Ta)
rho_sat = 0.622 / (Rair * Ta) * 610.78 * exp(21.875 * (Ta - Tm) / (Ta - 7.66))

! Assume that the rh varies according to a modification to 
! Pomeroy's humidity variation with height equation.
rh_offset = 1.0 - 0.027 * log(zRH)
sigma = (0.01 * RH - 1.0) * (rh_offset + 0.027 * log(z))
sigma = min(0.0,sigma)
sigma = max(-1.0,sigma)

alpha = 4.08 + 12.6 * z
rbar_r = 4.6e-5 * z**(-0.258)
xmbar = 4.0/3.0 * pi * rho_ice * rbar_r**3 * (1.0 + 3.0/alpha + 2.0/alpha**2)
rbar = ((3.0 * xmbar) / (4.0 * pi * rho_ice))**(0.33)
u_z = Utau/vkman * log(z/z0)
x_r = 0.005 * u_z**(1.36)
wbar = 1.1e7 * rbar**(1.8)

if (flag_susp) then

  ! Compute the sublimation loss rate coefficient for the suspension
  ! layer.
  V_r = wbar + 3.0 * x_r * cos(pi/4.0)
  xN_r = 2.0 * rbar * V_r / visc_air
  xNu = 1.79 + 0.606 * xN_r**(0.5)
  xSh = xNu
  tmp1 = (Ls * xM)/(Runi * Ta) - 1.0
  tmp2 = hcon_air * Ta * xNu
  top = 2.0 * pi * rbar * sigma
  bottom = Ls/tmp2 * tmp1 + 1.0/(D * rho_sat * xSh)
  V_susp = (top/bottom)/xmbar
  V_salt = 0.0

else

  ! Compute the sublimation loss rate coefficient for the saltation
  ! layer.
  V_rsalt = 0.68 * Utau + 2.3 * Utau_t
  xN_r = 2.0 * rbar * V_rsalt / visc_air
  xNu = 1.79 + 0.606 * xN_r**(0.5)
  xSh = xNu
  tmp1 = (Ls * xM)/(Runi * Ta) - 1.0
  tmp2 = hcon_air * Ta * xNu
  top = 2.0 * pi * rbar * sigma
  bottom = Ls/tmp2 * tmp1 + 1.0/(D * rho_sat * xSh)
  V_salt = (top/bottom)/xmbar
  V_susp = 0.0

end if

end subroutine SUBLIMATION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine COMPUTE_WIND_DIRECTION_ARRAYS(index_wdir_E,index_wdir_W,index_wdir_N,index_wdir_S,wind_WE,wind_SN)

use GRID, only: &
  Nx,Ny                   ! Grid dimensions

implicit none

integer, intent(inout) :: &
  index_wdir_E(Nx,2*Ny+1), &! Wind index array E
  index_wdir_W(Nx,2*Ny+1), &! Wind index array W
  index_wdir_N(Ny,2*Nx+1), &! Wind index array N
  index_wdir_S(Ny,2*Nx+1)   ! Wind index array S

real, intent(in) :: &
  wind_WE(Nx,Ny),          &! W->E (y) component of wind speed (m/s)
  wind_SN(Nx,Ny)            ! S->N (x) component of wind speed (m/s)

! Index whether the winds are blowing east or west.  The first
! column of the index array is the number of pairs of begining
! and ending array index of blocks of wind running in the same
! direction.

! Sweep looking for Westerly, Easterly, Southerly and Northerly winds
! wind_WE and wind_SN are passed in their natural (x, y) layout.
! swap_indices = .FALSE.: (row, position) = (x, y) -> W/E sweeps, transport along y
! swap_indices = .TRUE.:  (row, position) = (y, x) -> S/N sweeps, transport along x
call SWEEP_DIRECTION(index_wdir_W, wind_WE, Nx, Ny, .true.,  .false.)  ! Westerly:  wind_WE > 0
call SWEEP_DIRECTION(index_wdir_E, wind_WE, Nx, Ny, .false., .false.)  ! Easterly:  wind_WE < 0
call SWEEP_DIRECTION(index_wdir_S, wind_SN, Ny, Nx, .true.,  .true.)  ! Southerly: wind_SN > 0
call SWEEP_DIRECTION(index_wdir_N, wind_SN, Ny, Nx, .false., .true.)  ! Northerly: wind_SN < 0

end subroutine COMPUTE_WIND_DIRECTION_ARRAYS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SWEEP_DIRECTION(index_array, wind_data, dim1_size, dim2_size, &
                           positive_direction, swap_indices)

use GRID, only: &
  Nx,Ny                    ! Grid dimensions

implicit none

integer, intent(in) :: &
  dim1_size,               &! Size of the first dimension (rows)
  dim2_size                 ! Size of the second dimension (transport axis)

logical, intent(in) :: &
  positive_direction,      &! Boolean indicating if the wind direction is along increasing or decreasing index
  swap_indices              ! .FALSE.: W/E sweeps, (row, position) = (x, y)
                            ! .TRUE.:  S/N sweeps, (row, position) = (y, x)

integer, intent(out) :: &
  index_array(dim1_size, 2*dim2_size+1) ! Wind index array

real, intent(in) :: &
  wind_data(Nx,Ny)          ! Wind component in its natural (x, y) layout (m/s)

integer :: &
  dim1, dim2,              &! Row and transport-axis indexes
  npairs                    ! Number of pairs of block boundaries

real :: &
  wind_row(dim2_size),     &! Wind component along the transport axis for current row (m/s)
  sign_prev, sign_curr      ! Sign indicators of wind

do dim1 = 1, dim1_size

  ! Extract the wind component along the transport axis for this row,
  ! from the natural (x, y) layout of wind_data.
  do dim2 = 1, dim2_size
    if (swap_indices) then
      wind_row(dim2) = wind_data(dim2, dim1)   ! S/N: position = x (i), row = y (j)
    else
      wind_row(dim2) = wind_data(dim1, dim2)   ! W/E: row = x (i), position = y (j)
    end if
  end do

  ! Initialize first element sign
  sign_prev = sign(1.0, wind_row(1))

  ! Determine starting conditions based on target direction
  if ((positive_direction .and. sign_prev > 0.0) .or. &
      (.not. positive_direction .and. sign_prev < 0.0)) then
    npairs = 1
    index_array(dim1,2) = 1
  else
    npairs = 0
  end if

  ! Sweep through remaining elements in dim2 direction
  do dim2 = 2, dim2_size
    sign_curr = sign(1.0, wind_row(dim2))

    ! Check for sign change
    if (sign_curr /= sign_prev) then
      if ((positive_direction .and. sign_curr > 0.0) .or. &
          (.not. positive_direction .and. sign_curr < 0.0)) then
        ! Start of new target group
        npairs = npairs + 1
        index_array(dim1,npairs*2) = dim2
      else
        ! End of target group
        index_array(dim1,npairs*2+1) = dim2 - 1
      end if
    end if

    sign_prev = sign_curr
  end do

  ! Handle case where last element is still in target group
  sign_curr = sign(1.0, wind_row(dim2_size))
  if ((positive_direction .and. sign_curr > 0.0) .or. &
      (.not. positive_direction .and. sign_curr < 0.0)) then
    index_array(dim1,npairs*2+1) = dim2_size
  end if

  ! Store pair count
  index_array(dim1,1) = npairs
end do

end subroutine SWEEP_DIRECTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine UPDATE_MASS_CHANGES(Qs_WE,Qs_SN,dSWE_s, &
                               index_wdir_E,index_wdir_W,index_wdir_N,index_wdir_S, &
                               snowthickness,Ds_soft,snowdepth0,Sice0)

use GRID, only: &
  Nx,Ny                     ! Grid dimensions

implicit none

integer, intent(in) :: &
  index_wdir_E(Nx,2*Ny+1), &! Wind index array E
  index_wdir_W(Nx,2*Ny+1), &! Wind index array W
  index_wdir_N(Ny,2*Nx+1), &! Wind index array N
  index_wdir_S(Ny,2*Nx+1)   ! Wind index array S

real, intent(inout) :: &
  snowthickness(Nx,Ny),    &! Depth of the snowpack in the snow covered part, i.e. not scaled by fsnow (m)
  Qs_WE(Nx,Ny),             &! x component of either Qsalt or Qsusp (kg/m/s)
  Qs_SN(Nx,Ny),             &! y component of either Qsalt or Qsusp (kg/m/s)
  dSWE_s(Nx,Ny),           &! SWE change due to either saltation or suspension (kg/m^2)
  Ds_soft(Nx,Ny),          &! Soft snow thickness (m)
  snowdepth0(Nx,Ny),       &! Snow depth of snowdrift accumulation, averaged over the grid cell (m)
  Sice0(Nx,Ny)              ! Surface layer ice content for transported snow (kg/m^2)

! Sweep each wind direction

! Westerly winds (wnd_dir = 1) — j ascending, upwind = j-1, downwind = j+1
call COMPUTE_MASS_CHANGES_PER_DIRECTION(1, index_wdir_W, Nx, Ny, &
                                        dSWE_s, Qs_WE, snowthickness, Ds_soft, snowdepth0, Sice0)

! Easterly winds (wnd_dir = 2) — j descending, upwind = j+1, downwind = j-1
call COMPUTE_MASS_CHANGES_PER_DIRECTION(2, index_wdir_E, Nx, Ny, &
                                        dSWE_s, Qs_WE, snowthickness, Ds_soft, snowdepth0, Sice0)

! Southerly winds (wnd_dir = 3) — i ascending, upwind = i-1, downwind = i+1
call COMPUTE_MASS_CHANGES_PER_DIRECTION(3, index_wdir_S, Ny, Nx,  &
                                        dSWE_s, Qs_SN, snowthickness, Ds_soft, snowdepth0, Sice0)

! Northerly winds (wnd_dir = 4) — i descending, upwind = i+1, downwind = i-1
call COMPUTE_MASS_CHANGES_PER_DIRECTION(4, index_wdir_N, Ny, Nx, &
                                        dSWE_s, Qs_SN, snowthickness, Ds_soft, snowdepth0, Sice0)

end subroutine UPDATE_MASS_CHANGES

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine COMPUTE_MASS_CHANGES_PER_DIRECTION(wnd_dir, index_array, dim1_size, dim2_size, &
                                              dSWE_s, Qs_dir, snowthickness, Ds_soft, snowdepth0, Sice0)

use GRID, only: &
  Nx,Ny                                ! Grid dimensions

use DRIVING, only: &
  dt                                   ! Timestep (s)

use PARAMMAPS, only: &
  veg_shd                              ! Vegetation snow holding depth (m)

use STATE_VARIABLES, only: &
  Ds,                                 &! Snow layer thicknesses (m)
  fsnow,                              &! Snow cover fraction 
  Sice,                               &! Ice content of snow layers (kg/m^2)
  Sliq                                 ! Liquid content of snow layers (kg/m^2)

use PARAM_SNOWTRAN3D, only: &
  rho_snow                             ! Constant snowdrift deposition density (kg/m^3)

use LANDUSE, only: &
  cellsize,                           &! Grid cell size (m)
  tilefrac,                           &! Grid cell tile fraction - here = open fraction since SNOWTRAN3D only runs in open tile
  glacierfrac                          ! Glacier fraction

use FUNCTIONS_SNOWTRAN3D

implicit none

real, intent(inout) :: &
  dSWE_s(Nx,Ny),                      &! SWE change due to either saltation or suspension (kg/m^2)
  snowthickness(Nx,Ny),               &! Depth of the snowpack in the snow covered part, i.e. not scaled by fsnow (m)
  Ds_soft(Nx,Ny),                     &! Soft snow thickness (m)
  Qs_dir(Nx,Ny),                      &! x or y component of either Qsalt or Qsusp (kg/m/s)
  snowdepth0(Nx,Ny),                  &! Snow depth of snowdrift accumulation, averaged over the grid cell (m)
  Sice0(Nx,Ny)                         ! Surface layer ice content for transported snow (kg/m^2)

integer, intent(in) :: &
  wnd_dir,                            &! Wind direction index (1=West, 2=East, 3=South, 4=North)
  dim1_size,                          &! Size of the first dimension
  dim2_size                            ! Size of the second dimension

integer, intent(in) :: &
  index_array(dim1_size,2*dim2_size+1) ! Wind index array

real :: &
  dHS_gain,                          &! Snow depth gain due to either saltation or suspension (m)
  dHS_loss,                          &! Snow depth loss due to either saltation or suspension (m)
  dSWE_loc,                          &! SWE variation due to either saltation or suspension (kg/m^2)
  dSWE_loss,                         &! SWE loss due to either saltation or suspension (kg/m^2)
  Ds_hard,                           &! Hard snow thickness (m)
  snowdmin,                          &! Minimum snow depth to allow transport (m)
  swe_loc                             ! SWE in the snow covered part of the pixel (kg/m^2)

integer :: &
  kstart, kend,                      &! Index boundaries
  incr,                              &! Loop increment (1 or -1 depending on wind direction)
  dim1, dim2,                        &! Indexes of first and second dimensions
  k,                                 &! Temporary index
  ic, jc,                            &! x and y indexes of the considered grid cell
  iuw, juw,                           &! x and y indexes of the upwind grid cell
  idw, jdw                            ! x and y indexes of the downwind grid cell

real :: &
  openfrac(Nx,Ny)                      ! Grid cell open fraction (including glaciers)

! For snow transport, glaciers are considered open.
openfrac = tilefrac + glacierfrac

select case (wnd_dir)
case (1,3) ! Westerly and Southerly winds
  incr = 1
case (2,4) ! Easterly and Northerly winds
  incr = -1
end select

do dim1 = 1, dim1_size
  do k = 1, index_array(dim1,1)

    kstart = get_kstart(wnd_dir, dim1, k, index_array, dim1_size, dim2_size)
    kend   = get_kend(wnd_dir, dim1, k, index_array, dim1_size, dim2_size)

    do dim2 = kstart, kend, incr

      select case (wnd_dir)
      case (1,2) ! Westerly and Easterly winds
        ic = dim1
        jc = dim2
      case (3,4) ! Southerly and Northerly winds
        ic = dim2
        jc = dim1
      end select

      iuw = get_i_upwind(wnd_dir, ic)
      juw = get_j_upwind(wnd_dir, jc)
      idw = get_i_downwind(wnd_dir, ic)
      jdw = get_j_downwind(wnd_dir, jc)

      ! Stop transport at an internal wind reversal so the terminal cell deposits the incoming flux.
      if (dim2 == kend .and. idw >= 1 .and. idw <= Nx .and. &
          jdw >= 1 .and. jdw <= Ny) then
        Qs_dir(ic,jc) = 0.0
      end if

      ! Scaling by fsnow and openfrac is already done in SALTATION.
      dSWE_loc = dt * (Qs_dir(iuw,juw) - Qs_dir(ic,jc))  / cellsize(ic,jc)

      ! To conserve mass, we need to divide by the open fraction of the pixel
      ! to compensate for the final tile combination multiplying by this fraction -> / openfrac(ic,jc)
      if (openfrac(ic,jc) >= 0.1) then
        dSWE_loc = dSWE_loc / openfrac(ic,jc)
      else
        ! In theory, the flux is already zero but let's be safe.
        dSWE_loc = 0.0
      end if

      ! Case of a SWE loss
      if (dSWE_loc < epsilon(dSWE_loc)) then

        dSWE_loss = -dSWE_loc
        ! Make adjustments of snow loss for the case where there is no snow available
        ! on the ground (or captured within the vegetation) to be eroded.
        ! 1: Limit SWE loss to total SWE on the pixel
        swe_loc = sum(Sice(:,ic,jc) + Sliq(:,ic,jc))
        dSWE_loss = min(dSWE_loss,swe_loc)
        ! And adjust corresponding HS loss
        if (dSWE_loss > epsilon(dSWE_loss)) then
          call HS_FROM_SWE(dSWE_loss,dHS_loss,ic,jc)
        else
          dSWE_loss = 0.0
          dHS_loss = 0.0
        end if
        Ds_hard = snowthickness(ic,jc) - Ds_soft(ic,jc)
        snowdmin = max(veg_shd(ic,jc),Ds_hard)
        ! 2: Check if there is enough erodible HS (soft and out of vegetation)
        if (snowthickness(ic,jc) > snowdmin .and. fsnow(ic,jc) > epsilon(fsnow)) then
          ! 3: Check if the calculated HS loss is bigger than available erodible HS
          if (snowthickness(ic,jc) - dHS_loss / fsnow(ic,jc) <= snowdmin) then
            ! If so, limit HS loss to available erodible snow
            dHS_loss = (snowthickness(ic,jc) - snowdmin) * fsnow(ic,jc)
            ! And adjust corresponding SWE loss
            if (dHS_loss > epsilon(dHS_loss)) then
              call SWE_FROM_HS(dHS_loss,dSWE_loss,ic,jc)
            else
              dSWE_loss = 0.0
            end if
          end if
        else
          ! No available erodible snow: adjust SWE loss
          dSWE_loss = 0.0
        end if

        ! Proceed to actual ablation of eroded snow and update soft snow depth
        if (dSWE_loss > epsilon(dSWE_loss) .and. dHS_loss > epsilon(dHS_loss)) then
          call SNOW_ABLATION(dHS_loss,dSWE_loss,ic,jc)
          if (fsnow(ic,jc) > epsilon(fsnow)) then
            Ds_soft(ic,jc) = max(Ds_soft(ic,jc) - dHS_loss / fsnow(ic,jc),0.0)
          end if
          dSWE_loc = -dSWE_loss
          ! Finally, adjust outgoing flux of eroded snow
          Qs_dir(ic,jc) = Qs_dir(iuw,juw) + dSWE_loss * openfrac(ic,jc) * cellsize(ic,jc) / dt
        else
          dSWE_loc = 0.0
          Qs_dir(ic,jc) = Qs_dir(iuw,juw)
        end if

      ! Case of a SWE gain
      else

        ! Proceed to actual accumulation of transported snow
        dHS_gain = dSWE_loc / rho_snow
        if (dSWE_loc > epsilon(dSWE_loc) .and. dHS_gain > epsilon(dHS_gain)) then
          ! Add to the existing top layer.
          Sice0(ic,jc) = Sice0(ic,jc) + dSWE_loc
          snowdepth0(ic,jc) = snowdepth0(ic,jc) + dHS_gain
        else
          dSWE_loc = 0.0
          Qs_dir(ic,jc) = Qs_dir(iuw,juw)
        end if

      end if

      ! Update the net SWE change due to saltation or suspension transport for saving it later
      dSWE_s(ic,jc) = dSWE_s(ic,jc) + dSWE_loc

      ! Update the snow layer thicknesses
      snowthickness(ic,jc) = sum(Ds(:,ic,jc))

    end do
  end do
end do

end subroutine COMPUTE_MASS_CHANGES_PER_DIRECTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine COMPUTE_THRESHOLD_FRICTION_VELOCITY(Utau_t)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use PARAM_SNOWTRAN3D, only: &
  rho_snow            ! Constant snow density (kg/m^3)

use STATE_VARIABLES, only: &
  Nsnow,             &! Number of snow layers
  fsnow,             &! Snow cover fraction 
  Sice,              &! Ice content of snow layers (kg/m^2)
  Sliq,              &! Liquid content of snow layers (kg/m^2)
  Ds                  ! Snow layer thicknesses (m)

implicit none

real, intent(inout) :: &
  Utau_t(Nx,Ny)       ! Threshold friction velocity (m/s)

integer :: &
  i,j                 ! Point counters

real :: &
  rho_surf_snow       ! Density of top snow layer  (kg/m^3)


do j = 1, Ny
  do i = 1, Nx

    rho_surf_snow = rho_snow
    if (Nsnow(i,j) > 0 .and. Ds(1,i,j) > epsilon(Ds)) then
      rho_surf_snow = (Sice(1,i,j) + Sliq(1,i,j)) / Ds(1,i,j) / fsnow(i,j)
    end if

    ! Calculate the snow threshold friction velocity.
    if (rho_surf_snow <= 300.0) then
      Utau_t(i,j) = 0.10 * exp(0.003 * rho_surf_snow)
    else
      Utau_t(i,j) = 0.005 * exp(0.013 * rho_surf_snow)
    end if

  end do
end do

end subroutine COMPUTE_THRESHOLD_FRICTION_VELOCITY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
! Lateral redistribution of snow through gravity
! Implementation of the SnowSlide model (Bernhardt and Schulz, 2010)
! in the FSM model by Louis Quéno
! Reference: Quéno et al. (2024)
! ATTENTION: y is the W->E axis, while x is the S->N
! South of (i,j): (i-1,j)
! North of (i,j): (i+1,j)
! West of (i,j): (i,j-1)
! East of (i,j): (i,j+1)
!-----------------------------------------------------------------------
subroutine SNOWSLIDE(snowdepth0,Sice0,dSWE_slide)

use GRID, only: &
  Nx,Ny                     ! Grid dimensions

use STATE_VARIABLES, only: &
  fsnow,                   &! Snow cover fraction 
  Ds,                      &! Snow layer thicknesses (m)
  dSWE_tot_slide,          &! Cumulated SWE change due to snow slides (kg/m^2)
  index_sorted_dem          ! Location (i,j) of sorted grid points

use LANDUSE, only: &
  cellsize,                &! Grid cell size (m)
  dem,                     &! Terrain elevation (m)
  slope,                   &! Slope (deg)
  Shd,                     &! Snow holding depth (m)
  tilefrac,                &! Grid cell tile fraction - here = open fraction since SNOWSLIDE only runs in open tile
  glacierfrac               ! Glacier fraction

use PARAM_SNOWSLIDE, only: &
  dyn_ratio,               &! Dynamic snow holding depth ratio
  trig_ratio,              &! Hysteretic triggering snow holding depth ratio
  slope_min,               &! Minimum slope for snow slide occurrence (deg)
  Shd_min                   ! Minimum snow holding depth (m)

use FUNCTIONS_SNOWSLIDE

implicit none

real, intent(inout) :: &
  snowdepth0(Nx,Ny),       &! Snow depth of deposited snow (m)
  Sice0(Nx,Ny),            &! Ice content of deposited snow (kg/m^2)
  dSWE_slide(Nx,Ny)         ! SWE change due to snow slides (kg/m^2)

integer :: &
  i,j,                     &! x and y indexes of grid cell
  n,                       &! Vector point counter
  i_neigh,j_neigh,         &! x and y indexes of neighbouring grid cell
  k_neigh                   ! Neighbour index

logical :: &
  snow_excess_trigger,     &! Boolean to identify pixels where snow is available for avalanche triggering
  snow_excess_depo          ! Boolean to identify pixels where snow is available for avalanche continuation

logical :: &
  snow_depo(Nx,Ny)          ! Boolean to identify pixels where snow is deposited

real :: &
  delta,                   &! Distance between two neighbouring pixel centers (m)
  dswe,                    &! Mass of snow transported to a pixel (kg/m^2)
  elev,                    &! Elevation with snow depth of grid point n (m)
  elev_neigh,              &! Elevation with snow depth of neighboring grid point (m)
  weight_tot,              &! Total weight for normalization
  wt,                      &! Ratio of fresh avalanche deposit available for new slide
  rho_init,                &! Initial density of sliding snow (kg/m^3)
  rho_deposit,             &! Snow avalanche deposit density (kg/m^3)
  snowdepth_available,     &! Available snow depth for slide (m)
  snowdepth_available2,    &! Available snow depth for slide after substracting fresh deposit(m)
  swe_available,           &! Available SWE for slide (kg/m^2)
  swe_available2,          &! Available SWE for slide after substracting fresh deposit (kg/m^2)
  snowdepth_updated         ! Local snow depth updated in the snow slide loop (m)

real :: &
  Shd_corr(Nx,Ny),         &! Snow holding depth, thresholded (m)
  snowdepth(Nx,Ny),        &! Snow depth averaged over the grid cell before slides (m)
  openfrac(Nx,Ny)           ! Grid cell open fraction (including glaciers)

integer, parameter :: &
  N_neigh = 8               ! Number of neighbouring grid cells

integer :: &
  di(N_neigh),             &! Index offset in x for each neighbouring grid cell
  dj(N_neigh)               ! Index offset in y for each neighbouring grid cell

real :: &
  dfactor(N_neigh),        &! Distance factor for each neighbouring grid cell
  weight(N_neigh)           ! Weight of neighbouring grid cells for slide

! For snow transport, glaciers are considered open.
openfrac = tilefrac + glacierfrac

! Direction offset arrays (S, N, W, E, SW, SE, NW, NE)
di = (/ -1,  1,  0,  0, -1, -1,  1,  1 /)
dj = (/  0,  0, -1,  1, -1,  1, -1,  1 /)
dfactor = sqrt(real(di**2 + dj**2))

! Initialize snow_depo array
snow_depo(:,:) = .FALSE.

! Compute snowdepth before snow slides, and threshold the snow holding depth with Shd_min
do j = 1, Ny
  do i = 1, Nx
    snowdepth(i,j) = sum(Ds(:,i,j)) * fsnow(i,j)
    Shd_corr(i,j) = max(Shd(i,j),Shd_min)
  end do
end do

! Treating grid points from the highest to the lowest
do n = 1, Nx*Ny

  i = index_sorted_dem(n,1)
  j = index_sorted_dem(n,2)

  ! Start slide processes only if slope higher than the defined minimum. 
  ! If it is a pixel receiving avalanche snow, no slope threshold.
  if (slope(i,j) >= slope_min .or. snow_depo(i,j)) then

    ! Update snowdepth in case snow has been transported to this pixel earlier in the loop
    snowdepth_updated = sum(Ds(:,i,j)) * fsnow(i,j) + snowdepth0(i,j)

    ! Local elevation accounting for updated snowdepth
    elev = dem(i,j) + snowdepth_updated

    ! If an avalanche is occurring, the snow holding depth is reduced
    ! to mimic the dynamic effect
    if (snow_depo(i,j)) then
      Shd_corr(i,j) = Shd_corr(i,j) * dyn_ratio
      Shd_corr(i,j) = max(Shd_corr(i,j),Shd_min)
    end if

    ! Compute the depth of snow available for avalanche transport
    snowdepth_available = max(0.0, snowdepth_updated - Shd_corr(i,j))

    ! In absence of current avalanche, an hysteretic ratio of the snow holding depth must be overcome
    ! to trigger a new avalanche.
    snow_excess_trigger = (.not. snow_depo(i,j)) &
                          .and. (snowdepth_updated - Shd_corr(i,j) * trig_ratio > epsilon(snowdepth_available))
    ! If an avalanche is arriving at the pixel, the available snowdepth over holding threshold must be overcome 
    ! to continue the avalanche.
    snow_excess_depo = (snow_depo(i,j)) .and. (snowdepth_available > epsilon(snowdepth_available))

    ! There is a snow excess at the pixel.
    if (snow_excess_trigger .or. snow_excess_depo) then

      ! ---------------------------------------------------------------------------
      ! 1. Edge pixels. Excess snow is dumped out of the domain
      ! to avoid accumulation artefacts. It is not routed to domain pixels.
      ! ---------------------------------------------------------------------------
      if (i == 1 .or. i == Nx .or. j == 1 .or. j == Ny) then

        ! Move first snow coming from fresh avalanche deposit
        if (snowdepth0(i,j) - snowdepth_available > epsilon(snowdepth0)) then

          wt = snowdepth_available / snowdepth0(i,j)
          snowdepth0(i,j) = snowdepth0(i,j) - snowdepth_available
          swe_available = wt * Sice0(i,j)
          Sice0(i,j) = Sice0(i,j) - swe_available

        else if (snowdepth_available - snowdepth0(i,j) > epsilon(snowdepth0)) then

          snowdepth_available2 = snowdepth_available - snowdepth0(i,j)
          snowdepth0(i,j) = 0.0
          swe_available = Sice0(i,j)
          Sice0(i,j) = 0.0

          ! Compute the mass of snow available for transport
          call SWE_FROM_HS(snowdepth_available2,swe_available2,i,j)

          call SNOW_ABLATION(snowdepth_available2,swe_available2,i,j)

          swe_available = swe_available + swe_available2

        else ! snowdepth_available == snowdepth0(i,j)

          snowdepth0(i,j) = 0.0
          swe_available = Sice0(i,j)
          Sice0(i,j) = 0.0

        end if

        dSWE_slide(i,j) = dSWE_slide(i,j) - swe_available
        dSWE_tot_slide(i,j) = dSWE_tot_slide(i,j) - swe_available

      ! ---------------------------------------------------------------------------
      ! 2. Iner domain pixels. Excess snow is routed to lower pixels if there are.
      ! Mass transfers are weighted by local slopes. The forest fraction reduces
      ! the weight of the pixel. No snow is transferred to fully forested pixels.
      ! ---------------------------------------------------------------------------
      else

        ! Calculate weights for all 8 directions in a single loop        
        weight_tot = 0.0

        do k_neigh = 1, N_neigh

          i_neigh = i + di(k_neigh)
          j_neigh = j + dj(k_neigh)
          
          ! Elevation of neighbouring grid cell
          elev_neigh = dem(i_neigh,j_neigh) + sum(Ds(:,i_neigh,j_neigh)) * fsnow(i_neigh,j_neigh) + snowdepth0(i_neigh,j_neigh)
          ! Distance to neighbouring grid cell
          delta = (cellsize(i,j) + cellsize(i_neigh,j_neigh)) * dfactor(k_neigh) / 2.0
          ! Weight based on slope towards neighbour and open fraction of destination
          weight(k_neigh) = max(0.0, atan((elev - elev_neigh)/delta)) * openfrac(i_neigh,j_neigh)

          ! Accumulate total weight
          weight_tot = weight_tot + weight(k_neigh)

        end do
        
        ! if weight_tot == 0 then it's a "sink", no avalanche possible.
        if (weight_tot > epsilon(weight_tot)) then

          ! Normalize weights
          weight = weight / weight_tot

          ! Move first snow coming from fresh avalanche deposit
          if (snowdepth0(i,j) - snowdepth_available > epsilon(snowdepth0)) then

            wt = snowdepth_available / snowdepth0(i,j)
            snowdepth0(i,j) = snowdepth0(i,j) - snowdepth_available
            swe_available = wt * Sice0(i,j)
            Sice0(i,j) = Sice0(i,j) - swe_available

          else if (snowdepth_available - snowdepth0(i,j) > epsilon(snowdepth0)) then

            snowdepth_available2 = snowdepth_available - snowdepth0(i,j)
            snowdepth0(i,j) = 0.0
            swe_available = Sice0(i,j)
            Sice0(i,j) = 0.0

            ! Compute the mass of snow available for transport
            call SWE_FROM_HS(snowdepth_available2,swe_available2,i,j)

            call SNOW_ABLATION(snowdepth_available2,swe_available2,i,j)

            swe_available = swe_available + swe_available2

          else ! snowdepth_available == snowdepth0(i,j)

            snowdepth0(i,j) = 0.0
            swe_available = Sice0(i,j)
            Sice0(i,j) = 0.0

          end if

          ! Compute the density of the avalanche deposit depending on the initial snow density
          rho_init = swe_available / snowdepth_available ! snowdepth_available is > 0 anyway if snow_excess_trigger or snow_excess_depo
          rho_deposit = DEPOSIT_DENSITY(rho_init)

          dSWE_slide(i,j) = dSWE_slide(i,j) - swe_available
          dSWE_tot_slide(i,j) = dSWE_tot_slide(i,j) - swe_available

          ! Transport to neighbour pixels, weighted by elevation difference
          ! Higher pixels have a weight of 0
          do k_neigh = 1, N_neigh

            if (weight(k_neigh) > epsilon(weight(k_neigh))) then

              i_neigh = i + di(k_neigh)
              j_neigh = j + dj(k_neigh)
              
              dswe = weight(k_neigh) * swe_available
              ! The lower pixel only receives snow from the open part of the higher pixel.
              ! The SWE transfer is concentrated on the open part of the pixel, 
              ! hence division by open fraction to cancel counter-weighting later when combining tiles.
              dswe = dswe * openfrac(i,j) / openfrac(i_neigh,j_neigh)
              
              ! Update neighbour grid cell state
              Sice0(i_neigh,j_neigh) = Sice0(i_neigh,j_neigh) + dswe
              snowdepth0(i_neigh,j_neigh) = snowdepth0(i_neigh,j_neigh) + dswe / rho_deposit
              snow_depo(i_neigh,j_neigh) = .TRUE.
              dSWE_slide(i_neigh,j_neigh) = dSWE_slide(i_neigh,j_neigh) + dswe
              dSWE_tot_slide(i_neigh,j_neigh) = dSWE_tot_slide(i_neigh,j_neigh) + dswe

            endif

          end do

        end if

      end if

    end if

  end if

end do

end subroutine SNOWSLIDE

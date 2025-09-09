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
! Modified to be standalone - all variables passed as arguments
!-----------------------------------------------------------------------
subroutine SNOWSLIDE(Nx, Ny, Nsmax, snowdepth0, Sice0, dSWE_slide, &
                     fsnow, Ds, dSWE_tot_slide, index_sorted_dem, &
                     dem, slope, Shd, dyn_ratio, rho_deposit, &
                     slope_min, Shd_min, Ds_min, Tm, rho_snow, &
                     Sice, Sliq, Nsnow, Tsnow, histowet, &
                     rhos_min, rhos_max)

implicit none

! Input dimensions and parameters
integer, intent(in) :: Nx, Ny, Nsmax
real, intent(in) :: dyn_ratio, rho_deposit, slope_min, Shd_min
real, intent(in) :: Ds_min, Tm, rho_snow, rhos_min, rhos_max

! Input/output arrays
real, intent(inout) :: &
  snowdepth0(Nx,Ny),       &! Snow depth of deposited snow (m)
  Sice0(Nx,Ny),            &! Ice content of deposited snow (kg/m^2)
  dSWE_slide(Nx,Ny)         ! SWE change due to snow slides (kg/m^2)

! Input state arrays
real, intent(in) :: &
  fsnow(Nx,Ny),            &! Snow cover fraction
  Ds(Nsmax,Nx,Ny),         &! Snow layer thicknesses (m)
  dem(Nx,Ny),              &! Terrain elevation (m)
  slope(Nx,Ny),            &! Slope (deg)
  Shd(Nx,Ny)                ! Snow holding depth (m)

integer, intent(in) :: &
  index_sorted_dem(Nx*Ny,2), &! Location (i,j) of sorted grid points
  Nsnow(Nx,Ny)               ! Number of snow layers

real, intent(inout) :: &
  dSWE_tot_slide(Nx,Ny),   &! Cumulated SWE change due to snow slides (kg/m^2)
  Sice(Nsmax,Nx,Ny),       &! Ice content of snow layers (kg/m^2)
  Sliq(Nsmax,Nx,Ny),       &! Liquid content of snow layers (kg/m^2)
  Tsnow(Nsmax,Nx,Ny),      &! Snow layer temperatures (K)
  histowet(Nsmax,Nx,Ny)     ! Historical variable for past wetting of a layer (0-1)

integer :: &
  i,j,                     &! Point counters
  n                         ! Vector point counter

logical :: &
  snow_depo(Nx,Ny)          ! Boolean to identify pixels where snow is deposited

real :: &
  dswe,                    &! Mass of snow transported to a pixel (kg/m2)
  elev,                    &! Elevation with snowdepth of grid point n (m)
  elev_S,                  &! Elevation with snow depth of South grid point (m)
  elev_N,                  &! Elevation with snow depth of North grid point (m)
  elev_W,                  &! Elevation with snow depth of West grid point (m)
  elev_E,                  &! Elevation with snow depth of East grid point (m)
  elev_SW,                 &! Elevation with snow depth of South-West grid point (m)
  elev_SE,                 &! Elevation with snow depth of South-East grid point (m)
  elev_NW,                 &! Elevation with snow depth of North-West grid point (m)
  elev_NE,                 &! Elevation with snow depth of North-East grid point (m)
  w_S,                     &! Weight of South grid point for slide
  w_N,                     &! Weight of North grid point for slide
  w_W,                     &! Weight of West grid point for slide
  w_E,                     &! Weight of East grid point for slide
  w_SW,                    &! Weight of South-West grid point for slide
  w_SE,                    &! Weight of South-East grid point for slide
  w_NW,                    &! Weight of North-West grid point for slide
  w_NE,                    &! Weight of North-East grid point for slide
  wt,                      &! Ratio of fresh avalanche deposit available for new slide
  delev_tot,               &! Total elevation difference with lower pixels (m)
  snowdepth_available,     &! Available snow depth for slide (m)
  snowdepth_available2,    &! Available snow depth for slide after substracting fresh deposit(m)
  swe_available,           &! Available SWE for slide (kg/m2)
  swe_available2,          &! Available SWE for slide after substracting fresh deposit (kg/m2)
  snowdepth_updated         ! Local snow depth updated in the snow slide loop (m)

real :: &
  Shd_corr(Nx,Ny),         &! Snow holding depth, thresholded (m)
  snowdepth(Nx,Ny)          ! Snow depth averaged over the grid cell before slides (m)

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
    ! to micic the dynamic effect
    if (snow_depo(i,j)) then
      Shd_corr(i,j) = Shd_corr(i,j) * dyn_ratio
      Shd_corr(i,j) = max(Shd_corr(i,j),Shd_min)
    end if

    ! Compute the depth of snow available for avalanche transport
    snowdepth_available = max(0.0, snowdepth_updated - Shd_corr(i,j))

    ! if snowdepth_available == 0 then there is no snow available to slide.

    ! There is a snow excess at the pixel.
    if (snowdepth_available > epsilon(snowdepth_available)) then

      if (i == 1 .or. i == Nx .or. j == 1 .or. j == Ny) then
        ! Case 1: edge pixels. Excess snow is dumped out of the domain
        ! to avoid accumulation artefacts. It is not routed to domain pixels.

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
          call SWE_FROM_HS(snowdepth_available2, swe_available2, i, j, &
                           Nsmax, Nx, Ny, Nsnow, fsnow, Sice, Sliq, Ds, &
                           rhos_min, rhos_max, rho_snow)

          call SNOW_ABLATION(snowdepth_available2, swe_available2, i, j, &
                             Nsmax, Nx, Ny, Sice, Sliq, Ds, histowet, &
                             Nsnow, fsnow, Tsnow, Ds_min, Tm)

          swe_available = swe_available + swe_available2

        else ! snowdepth_available == snowdepth0(i,j)

          snowdepth0(i,j) = 0.0
          swe_available = Sice0(i,j)
          Sice0(i,j) = 0.0

        end if

        dSWE_slide(i,j) = dSWE_slide(i,j) - swe_available
        dSWE_tot_slide(i,j) = dSWE_tot_slide(i,j) - swe_available

      else
        ! Case 2: inner domain pixels. Excess snow is routed to lower pixels if there are.
        ! Mass transfers are weighted by elevation (snowdepth included) differences.

        ! South of (i,j): (i-1,j)
        elev_S = dem(i-1,j) + sum(Ds(:,i-1,j)) * fsnow(i-1,j) + snowdepth0(i-1,j)
        w_S = max(0.0, elev - elev_S)
        ! North of (i,j): (i+1,j)
        elev_N = dem(i+1,j) + sum(Ds(:,i+1,j)) * fsnow(i+1,j) + snowdepth0(i+1,j)
        w_N = max(0.0, elev - elev_N)
        ! West of (i,j): (i,j-1)
        elev_W = dem(i,j-1) + sum(Ds(:,i,j-1)) * fsnow(i,j-1) + snowdepth0(i,j-1)
        w_W = max(0.0, elev - elev_W)
        ! East of (i,j): (i,j+1)
        elev_E = dem(i,j+1) + sum(Ds(:,i,j+1)) * fsnow(i,j+1) + snowdepth0(i,j+1)
        w_E = max(0.0, elev - elev_E)
        ! South-West of (i,j): (i-1,j-1)
        elev_SW = dem(i-1,j-1) + sum(Ds(:,i-1,j-1)) * fsnow(i-1,j-1) + snowdepth0(i-1,j-1)
        w_SW = max(0.0, elev - elev_SW)
        ! South-East of (i,j): (i-1,j+1)
        elev_SE = dem(i-1,j+1) + sum(Ds(:,i-1,j+1)) * fsnow(i-1,j+1) + snowdepth0(i-1,j+1)
        w_SE = max(0.0, elev - elev_SE)
        ! North-West of (i,j): (i+1,j-1)
        elev_NW = dem(i+1,j-1) + sum(Ds(:,i+1,j-1)) * fsnow(i+1,j-1) + snowdepth0(i+1,j-1)
        w_NW = max(0.0, elev - elev_NW)
        ! North-East of (i,j): (i+1,j+1)
        elev_NE = dem(i+1,j+1) + sum(Ds(:,i+1,j+1)) * fsnow(i+1,j+1) + snowdepth0(i+1,j+1)
        w_NE = max(0.0, elev - elev_NE)

        delev_tot = w_S + w_N + w_W + w_E + w_SW + w_SE + w_NW + w_NE
        ! if delev_tot == 0 then it's a "sink", no avalanche possible.

        if (delev_tot > epsilon(delev_tot)) then

          w_S = w_S / delev_tot
          w_N = w_N / delev_tot
          w_W = w_W / delev_tot
          w_E = w_E / delev_tot
          w_SW = w_SW / delev_tot
          w_SE = w_SE / delev_tot
          w_NW = w_NW / delev_tot
          w_NE = w_NE / delev_tot

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
            call SWE_FROM_HS(snowdepth_available2, swe_available2, i, j, &
                             Nsmax, Nx, Ny, Nsnow, fsnow, Sice, Sliq, Ds, &
                             rhos_min, rhos_max, rho_snow)

            call SNOW_ABLATION(snowdepth_available2, swe_available2, i, j, &
                               Nsmax, Nx, Ny, Sice, Sliq, Ds, histowet, &
                               Nsnow, fsnow, Tsnow, Ds_min, Tm)

            swe_available = swe_available + swe_available2

          else ! snowdepth_available == snowdepth0(i,j)

            snowdepth0(i,j) = 0.0
            swe_available = Sice0(i,j)
            Sice0(i,j) = 0.0

          end if

          dSWE_slide(i,j) = dSWE_slide(i,j) - swe_available
          dSWE_tot_slide(i,j) = dSWE_tot_slide(i,j) - swe_available

          ! Transport to neighbour pixels, weighted by elevation difference
          ! Higher pixels have a weight of 0

          ! South
          if (w_S > epsilon(w_S)) then
            dswe = w_S * swe_available
            Sice0(i-1,j) = Sice0(i-1,j) + dswe
            snowdepth0(i-1,j) = snowdepth0(i-1,j) + dswe / rho_deposit
            snow_depo(i-1,j) = .TRUE.
            dSWE_slide(i-1,j) = dSWE_slide(i-1,j) + dswe
            dSWE_tot_slide(i-1,j) = dSWE_tot_slide(i-1,j) + dswe
          end if

          ! North
          if (w_N > epsilon(w_N)) then
            dswe = w_N * swe_available
            Sice0(i+1,j) = Sice0(i+1,j) + dswe
            snowdepth0(i+1,j) = snowdepth0(i+1,j) + dswe / rho_deposit
            snow_depo(i+1,j) = .TRUE.
            dSWE_slide(i+1,j) = dSWE_slide(i+1,j) + dswe
            dSWE_tot_slide(i+1,j) = dSWE_tot_slide(i+1,j) + dswe
          end if

          ! West
          if (w_W > epsilon(w_W)) then
            dswe = w_W * swe_available
            Sice0(i,j-1) = Sice0(i,j-1) + dswe
            snowdepth0(i,j-1) = snowdepth0(i,j-1) + dswe / rho_deposit
            snow_depo(i,j-1) = .TRUE.
            dSWE_slide(i,j-1) = dSWE_slide(i,j-1) + dswe
            dSWE_tot_slide(i,j-1) = dSWE_tot_slide(i,j-1) + dswe
          end if

          ! East
          if (w_E > epsilon(w_E)) then
            dswe = w_E * swe_available
            Sice0(i,j+1) = Sice0(i,j+1) + dswe
            snowdepth0(i,j+1) = snowdepth0(i,j+1) + dswe / rho_deposit
            snow_depo(i,j+1) = .TRUE.
            dSWE_slide(i,j+1) = dSWE_slide(i,j+1) + dswe
            dSWE_tot_slide(i,j+1) = dSWE_tot_slide(i,j+1) + dswe
          end if

          ! South-West
          if (w_SW > epsilon(w_SW)) then
            dswe = w_SW * swe_available
            Sice0(i-1,j-1) = Sice0(i-1,j-1) + dswe
            snowdepth0(i-1,j-1) = snowdepth0(i-1,j-1) + dswe / rho_deposit
            snow_depo(i-1,j-1) = .TRUE.
            dSWE_slide(i-1,j-1) = dSWE_slide(i-1,j-1) + dswe
            dSWE_tot_slide(i-1,j-1) = dSWE_tot_slide(i-1,j-1) + dswe
          end if

          ! South-East
          if (w_SE > epsilon(w_SE)) then
            dswe = w_SE * swe_available
            Sice0(i-1,j+1) = Sice0(i-1,j+1) + dswe
            snowdepth0(i-1,j+1) = snowdepth0(i-1,j+1) + dswe / rho_deposit
            snow_depo(i-1,j+1) = .TRUE.
            dSWE_slide(i-1,j+1) = dSWE_slide(i-1,j+1) + dswe
            dSWE_tot_slide(i-1,j+1) = dSWE_tot_slide(i-1,j+1) + dswe
          end if

          ! North-West
          if (w_NW > epsilon(w_NW)) then
            dswe = w_NW * swe_available
            Sice0(i+1,j-1) = Sice0(i+1,j-1) + dswe
            snowdepth0(i+1,j-1) = snowdepth0(i+1,j-1) + dswe / rho_deposit
            snow_depo(i+1,j-1) = .TRUE.
            dSWE_slide(i+1,j-1) = dSWE_slide(i+1,j-1) + dswe
            dSWE_tot_slide(i+1,j-1) = dSWE_tot_slide(i+1,j-1) + dswe
          end if

          ! North-East
          if (w_NE > epsilon(w_NE)) then
            dswe = w_NE * swe_available
            Sice0(i+1,j+1) = Sice0(i+1,j+1) + dswe
            snowdepth0(i+1,j+1) = snowdepth0(i+1,j+1) + dswe / rho_deposit
            snow_depo(i+1,j+1) = .TRUE.
            dSWE_slide(i+1,j+1) = dSWE_slide(i+1,j+1) + dswe
            dSWE_tot_slide(i+1,j+1) = dSWE_tot_slide(i+1,j+1) + dswe
          end if

        end if

      end if

    end if

  end if

end do

end subroutine SNOWSLIDE

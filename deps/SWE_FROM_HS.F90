!-----------------------------------------------------------------------
! Equivalent SWE of surface HS
! Modified to be standalone - all variables passed as arguments
!-----------------------------------------------------------------------
subroutine SWE_FROM_HS(hs, swe, i, j, Nsmax, Nx, Ny, &
                       Nsnow, fsnow, Sice, Sliq, Ds, &
                       rhos_min, rhos_max, rho_snow)

! This subroutine calculate the average SWE of a given snow depth at the top of the snowpack (e.g. to be eroded)

implicit none

! Input parameters and dimensions
integer, intent(in) :: Nsmax, Nx, Ny, i, j
real, intent(in) :: hs, rhos_min, rhos_max, rho_snow

! Input state arrays
integer, intent(in) :: &
  Nsnow(Nx,Ny)          ! Number of snow layers
real, intent(in) :: &
  fsnow(Nx,Ny),        &! Snow cover fraction
  Sice(Nsmax,Nx,Ny),   &! Ice content of snow layers (kg/m^2)
  Sliq(Nsmax,Nx,Ny),   &! Liquid content of snow layers (kg/m^2)
  Ds(Nsmax,Nx,Ny)       ! Snow layer thicknesses (m)

! Output
real, intent(out) :: &
  swe                   ! SWE on the top of the snowpack (kg/m^2)

real :: &
  dDs,               &! Snow thickness on the top of the snowpack (m)
  rho_avg,           &! Bulk snow density (kg/m^3)
  weight,            &! Layer weight
  Ds_tmp,            &! Temporary snow thickness variable (m)
  eps,               &! Epsilon Ds tolerance to avoid instabilities
  swe_tmp,           &! Temporary SWE variable (kg/m^2)
  swe_layer,         &! Layer SWE (kg/m^2)
  snowthickness       ! Depth of the snowpack in the snow covered part, i.e. not scaled by fsnow (m)

integer :: &
  k                   ! Snow layer counters

! Initialization
eps = 1e-4
rho_avg = rho_snow
swe = 0.0

if (fsnow(i,j) > epsilon(fsnow)) then
  
  dDs = hs / fsnow(i,j)
  snowthickness = sum(Ds(:,i,j))

  if (dDs > epsilon(dDs) .and. dDs <= snowthickness + eps) then

    if (dDs >= snowthickness) then
    
      ! Fix computing approximations to avoid instabilities
      swe = sum(Sice(:,i,j)+Sliq(:,i,j))
      
    else
    
      ! Normal case   
      k = 1
      swe_tmp = 0.0
      Ds_tmp = 0.0

      do while (k <= Nsnow(i,j) .and. Ds_tmp < dDs)
        swe_layer = Sice(k,i,j) + Sliq(k,i,j)
        if (dDs - Ds_tmp > Ds(k,i,j)) then
          Ds_tmp = Ds_tmp + Ds(k,i,j)
          swe_tmp = swe_tmp + swe_layer
        else
          weight = (dDs-Ds_tmp)/Ds(k,i,j)
          swe_tmp = swe_tmp + weight*swe_layer
          Ds_tmp = dDs
        end if
        k = k+1
      end do

      swe = swe_tmp

    end if

    rho_avg = swe / hs

  else if (dDs < -epsilon(dDs)) then

    write(*,*) 'WARNING SWE_FROM_HS: dDs < 0. CHECK/DEBUG!', i, j
    write(*,*) 'dDs', dDs
    error stop
    
  else if (dDs > snowthickness + eps) then

    write(*,*) 'WARNING SWE_FROM_HS: dDs > snowthickness. CHECK/DEBUG!', i, j
    write(*,*) 'dDs', dDs, 'snowthickness', snowthickness, 'dDs-snowthickness', dDs-snowthickness
    error stop

  end if

end if

if ((rho_avg < rhos_min - 0.5 .or. rho_avg > rhos_max + 0.5) .and. Ds_tmp > 0.001) then
  write(*,*) 'WARNING SWE_FROM_HS: invalid density.', rho_avg, i, j
  write(*,*) 'Ds: ', Ds(:,i,j)
  write(*,*) 'SWE: ', Sice(:,i,j) + Sliq(:,i,j)
  write(*,*) 'rho: ', (Sice(1:Nsnow(i,j),i,j) + Sliq(1:Nsnow(i,j),i,j)) / Ds(1:Nsnow(i,j),i,j) / fsnow(i,j)
  error stop
end if

end subroutine SWE_FROM_HS

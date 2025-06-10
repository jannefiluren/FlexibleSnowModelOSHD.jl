!-----------------------------------------------------------------------
! Surface and canopy net shortwave radiation
!-----------------------------------------------------------------------
subroutine RADIATION(alb,SWsrf,SWveg,Sdirt,Sdift,asrf_out,SWsci,LWt)

use MODCONF, only: CANMOD, RADSBG, ALBEDO, OSHDTN,ALRADT

use MODPERT, only: ALPERT

use MODTILE, only: tthresh 

use CONSTANTS, only: &
  I0,                &! Solar constant (W/m^2)
  sb,                &! Stefan-Boltzmann constant (W/m^2/K^4)
  Tm                  ! Melting point (K)

use DRIVING, only: &
  year,              &! Year
  month,             &! Month of year
  day,               &! Day of month
  hour,              &! Hour of day
  dt,                &! Timestep (s)
  LW,                &! Incoming longwave radiation (W/m2)
  Sdif,              &! Diffuse shortwave radiation (W/m^2)
  Sdir,              &! Direct-beam shortwave radiation (W/m^2)
  Sf,                &! Snowfall rate (kg/m2/s)
  Sf24h,             &! Snowfall 24hr (kg/m2)
  Ta,                &! Air temperature (K)
  Tv,                &! Time-varying transmissivity for dSWR
  alP,               &! Albedo perturbation for fresh snow
  Sdird               ! Direct-beam shortwave radiation, per horizontal surface area (W/m2)

use GRID, only: &
  Nx,Ny               ! Grid dimensions

use PARAMETERS, only: &
  adm,               &! Melting snow albedo decay time (h)
  asmx,              &! Maximum albedo for fresh snow
  asmn,              &! Minimum albedo for melting snow
  avg0,              &! Snow-free vegetation albedo
  avgs,              &! Snow-covered vegetation albedo
  Talb,              &! Albedo decay temperature threshold (C)
  tcld,              &! Cold snow albedo decay timescale (s)
  tmlt,              &! Melting snow albedo decay timescale (s)
  adfs,              &! Forest albedo decay factor
  adfl,              &! Forest albedo decay factor
  fsar,              &! Albedo adjustment range for canopy dependence 
  Sfmin               ! Minimum 24h snowfall to refresh albedo (kg/m^2)

use PARAMMAPS, only: &
  adc,               &! Cold snow albedo decay time (h)
  afs,               &! Maximum albedo for fresh snow
  alb0,              &! Snow-free ground albedo
  fsky,              &! Sky view fraction
  scap,              &! Canopy snow capacity (kg/m^2)
  trcn                ! Canopy transmissivity

use STATE_VARIABLES, only: &
  albs,              &! Snow albedo
  Sice,              &! Ice content of snow layers (kg/m^2)
  Sliq,              &! Liquid content of snow layers (kg/m^2)
  fsnow,             &! Snow cover fraction
  Sveg,              &! Snow mass on vegetation (kg/m^2)
  Tsrf                ! Surface skin temperature (K)

use LANDUSE, only: &
  fsky_terr,         &! Terain sky view factor
  fveg,              &! Canopy cover fraction
  tilefrac            ! Grid cell tile fraction

implicit none

real, intent(inout) :: &
  alb(Nx,Ny),        &! Albedo
  Sdirt(Nx,Ny),      &! Incoming direct beam radiation corrected for subgrid topography (W/m^2)
  Sdift(Nx,Ny),      &! Incoming diffuse radiation corrected for subgrid topography (W/m^2)
  SWsrf(Nx,Ny),      &! Net SW radiation absorbed by the surface (W/m^2)
  SWveg(Nx,Ny)        ! Net SW radiation absorbed by vegetation (W/m^2)

integer :: &
  i,j                 ! Point counters

real, intent(out) :: &
  asrf_out(Nx,Ny),   &! Surface albedo
  LWt(Nx,Ny),        &! Incoming longwave radiation corrected for subgrid topography (W/m^2)
  SWsci(Nx,Ny)        ! Subcanopy incoming SWR (W/m^2)

real :: &
  adc_loc,           &! Local cold snow albedo decay time (h)
  adm_loc,           &! Local melting snow albedo decay time (h)
  afs_loc,           &! Local maximum albedo for fresh snow
  alim,              &! Limiting snow albedo
  acan,              &! Canopy albedo
  asrf,              &! Surface albedo
  aveg,              &! Vegetation albedo
  fcans,             &! Canopy snowcover fraction
  rt,                &! Reciprocal timescale for albedo adjustment (1/s)
  SWEtmp,            &! Temporary snow water equivalent for snow covered fraction calculation (kg/m^2)
  Sdif_aux,          &! Auxiliary variable containing diffuse SW radiation relevant to canopy transfer
  tau,               &! Snow albedo decay timescale (s)
  tdif,              &! Canopy transmissivity for diffuse radiation
  tdir,              &! Canopy transmissivity for direct-beam radiation
  SWtopo_out,        &! Outgoing SW radiation corrected for subgrid topography (W/m^2)
  Sun_elev            ! Solar elevation angle

! Snow albedo
do j = 1, Ny
do i = 1, Nx

  if (tilefrac(i,j) < tthresh) goto 1 ! exclude points outside tile of interest
  
  ! LQ we need to copy `adc(i,j)` et al. into `*_loc` to avoid erroneously updating it every  timestep.
  adc_loc = adc(i,j)
  adm_loc = adm
  afs_loc = afs(i,j)

  if (ALBEDO == 0) then
  ! Diagnostic
    albs(i,j) = asmn + (afs_loc - asmn)*(Tsrf(i,j) - Tm) / Talb
    if (albs(i,j) < min(afs_loc, asmn)) albs(i,j) = min(afs_loc, asmn)
    if (albs(i,j) > max(afs_loc, asmn)) albs(i,j) = max(afs_loc, asmn)

  elseif (ALBEDO == 1) then
  ! Prognostic
    tau = tcld
    if (Tsrf(i,j) >= Tm) tau = tmlt
    ! Forest adjustments -> not yet properly tested for OSHD but option currently unused
    if (month > 4 .AND. month < 10) then
      tau = 70*3600
    end if
    if (fveg(i,j) > 0 .and. Sdir(i,j) > epsilon(Sdir(i,j))) then
      tau = tau /((1-trcn(i,j)*fsky(i,j))*(1+adfl*Tv(i,j)) + adfs*Tv(i,j))
    else if (fveg(i,j) > 0 .and. Sdif(i,j) > epsilon(Sdif(i,j))) then 
      tau = tau/((1-trcn(i,j)*fsky(i,j)) + adfs*trcn(i,j)*fsky(i,j))
    else if (fveg(i,j) > 0 .and. (Sdir(i,j) + Sdif(i,j) <= epsilon(Sdir(i,j)+Sdif(i,j)))) then
      tau = tau/(2-trcn(i,j)*fsky(i,j))
    end if     

    rt = 1/tau + Sf(i,j)/Sfmin
    alim = (asmn/tau + Sf(i,j)*afs_loc/Sfmin)/rt
    albs(i,j) = alim + (albs(i,j) - alim)*exp(-rt*dt)
    if (albs(i,j) < min(afs_loc, asmn)) albs(i,j) = min(afs_loc, asmn)
    if (albs(i,j) > max(afs_loc, asmn)) albs(i,j) = max(afs_loc, asmn)

  else ! ALBEDO == 2
  ! Prognostic, tuned, copied from JIM
    SWEtmp = sum(Sice(:,i,j) + Sliq(:,i,j))

    ! BC 08.23: aspect-dependent albedo tuning. Activated for oper season 2024
    ! or optionally.
    if ((ALRADT == 1) .OR. (OSHDTN == 1)) then
      ! BC Oct 23: Jan's suggestion: modify only when the decay rate should be increased (ad* DECREASE), not decreased
      if ((Sdir(i,j) > epsilon(Sdir(i,j))) .AND. (Sdird(i,j) < Sdir(i,j))) then
        adm_loc = adm_loc * (Sdird(i,j))/(Sdir(i,j))
        adc_loc = adc_loc * (Sdird(i,j))/(Sdir(i,j))
        if (adm_loc < epsilon(adm_loc)) then
          adm_loc = epsilon(adm_loc)
        endif
        if (adc_loc < epsilon(adc_loc)) then
          adc_loc = epsilon(adc_loc)
        endif
      endif
    endif

    if (ALPERT) then
      adm_loc = adm_loc * alP(i,j)
      adc_loc = adc_loc * alP(i,j)
    endif
    if (Tsrf(i,j) >= Tm) then
      albs(i,j) = (albs(i,j) - asmn)*exp(-(dt/3600)/adm_loc) + asmn
    else
      albs(i,j) = albs(i,j) - (dt/3600)/adc_loc
    end if
    if (SWEtmp < 75.0) then ! more stuff showing on and up through snow
      afs_loc = afs_loc * 0.80
    end if
    ! Reset to fresh snow albedo (wasn't originally available; only else term)
    if ((Sf(i,j) * dt) > 0.0 .AND. Sf24h(i,j) > Sfmin) then
      albs(i,j) = afs_loc
    else
      albs(i,j) = albs(i,j) + (afs_loc - albs(i,j))*Sf(i,j)*dt/Sfmin
    end if
    !! End Adjustments
    if (albs(i,j) > afs_loc) albs(i,j) = afs_loc
    if (albs(i,j) < asmn) albs(i,j) = asmn
  endif

  1 continue ! Exclude points  
  
end do
end do

! Surface and canopy net shortwave radiation
do j = 1, Ny
do i = 1, Nx

  if (tilefrac(i,j) < tthresh) goto 2 ! exclude points outside tile of interest

  ! Surface albedo
  asrf = albs(i,j)*(1-fveg(i,j)*fsar)  
  if (fsnow(i,j) <= epsilon(fsnow)) then
    asrf = alb0(i,j)
    albs(i,j) = alb0(i,j)
  endif

  ! Partial snowcover on canopy
  fcans = 0
  if (scap(i,j) > epsilon(scap)) fcans = Sveg(i,j) / scap(i,j)
  aveg = (1 - fcans)*avg0 + fcans*avgs
  acan = fveg(i,j)*aveg
  ! Canopy surface albedo for computing terrain radiation over canopy
  alb(i,j) = fveg(i,j)*aveg + (1-fveg(i,j))*asrf
  
  ! Surface albedo is stored in asurf_out to write in results
  asrf_out(i,j) = alb(i,j)
  
  if (RADSBG == 1) then
    ! Call Subgrid parameterization for SW radiation to compute SWtopo,netto SWtn
    call SWRADTOPO(alb(i,j),Sdir(i,j),Sdif(i,j),SWsrf(i,j),Sdirt(i,j),Sdift(i,j),SWtopo_out,Sun_elev,year,month,day,hour,i,j)
  endif

  if (RADSBG == 0) then
    SWtopo_out =  alb(i,j)*(Sdir(i,j)+Sdif(i,j))
    Sdirt(i,j) = Sdir(i,j)
    Sdift(i,j) = Sdif(i,j) 
  endif

  ! Solar radiation trasmission 
  if (CANMOD == 0) then
    SWveg(i,j) = 0
    SWsrf(i,j) = (1 - alb(i,j))*(Sdir(i,j)+Sdif(i,j))
    SWsci(i,j) = Sdift(i,j)+Sdir(i,j)  
  endif

  if (CANMOD == 1) then
    Sdif_aux = fsky(i,j)/fsky_terr(i,j)*Sdift(i,j)
    tdif = trcn(i,j)
    tdir = Tv(i,j) 

    ! Effective albedo and net radiation
    alb(i,j) = acan + (1 - acan)*asrf*tdif**2
    if (Sdif_aux + Sdirt(i,j) > epsilon(Sdift(i,j)))  & 
      alb(i,j) = (acan*(Sdif_aux+tdir*Sdirt(i,j)) + asrf*tdif*(tdif*Sdif_aux+tdir*Sdirt(i,j))) / &
                (Sdif_aux + Sdirt(i,j))
    SWsrf(i,j) = (1 - asrf)*(tdif*Sdif_aux + tdir*Sdirt(i,j))
    SWveg(i,j) = ((1-tdif)*(1-aveg)+tdif*asrf*(1-tdif))*Sdif_aux + &   
                  (tdir*fveg(i,j)*(1-aveg)+tdir*asrf*(1-tdif))*Sdirt(i,j)   ! local SWR absorption by vegetation correlates with local tdir  
    SWsci(i,j) = tdif*Sdif_aux + tdir*Sdirt(i,j)
  endif 
  
  ! Thermal emissions from surroundings
  ! Terrain LWR if not calculated later
  LWt(i,j) = fsky_terr(i,j)*LW(i,j) + (1 - fsky_terr(i,j))*sb*Ta(i,j)**4

  ! LW overwritten by LWt only if EBALFOR is not used, where terrain impacts are accounted for already
  if (CANMOD == 0 .OR. fveg(i,j) == 0) then
    LW(i,j) = LWt(i,j)
  endif
            
  2 continue 
  
end do
end do

end subroutine RADIATION

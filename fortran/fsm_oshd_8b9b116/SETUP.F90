!-----------------------------------------------------------------------
! Set parameters, initialize prognostic variables and write metadata
!-----------------------------------------------------------------------
subroutine SETUP

use MODCONF, only: ALBEDO, CANMOD, CONDCT, DENSTY, EXCHNG, HYDROL, &
                   SNFRAC, RADSBG, ZOFFST, OSHDTN, HN_ON, FOR_HN, ALRADT, &
                   SNTRAN, SNSLID, SNOLAY, CHECKS

use MODPERT, only: Z0PERT, WCPERT, FSPERT, ALPERT, SLPERT

use MODTILE, only: TILE, tthresh 

use MODOUTPUT, only: LIST_DIAG_RESULTS, LIST_STATE_RESULTS

use CONSTANTS

use DIAGNOSTICS, only : Nave

use DRIVING, only: &
  dt,                &! Timestep (s)
  LW,                &! Incoming longwave radiation (W/m2)
  Ps,                &! Surface pressure (Pa)
  Qa,                &! Specific humidity (kg/kg)
  RH,                &! Relative humidity (%)
  Rf,                &! Rainfall rate (kg/m2/s)
  Sf,                &! Snowfall rate (kg/m2/s)
  Sf24h,             &! Snowfall 24hr (kg/m2)
  Sdir,              &! Incoming direct beam radiation on flat,unobstructed surface (W/m2)
  Sdif,              &! Incoming diffuse radiation on flat,unobstructed (W/m2)
  Ta,                &! Air temperature (K)
  Tv,                &! Time-varying canopy transmissivity for dSWR (-)
  Ua,                &! Wind speed (m/s)
  Udir,              &! Wind direction (degrees, clockwise from N)
  zT,                &! Temperature measurement height (m)
  zU,                &! Wind speed measurement height (m)
  zRH,               &! Relative humidity measurement height (m)
  z0P,               &! z0 perturbations
  wcP,               &! liquid water capacity perturbations
  fsP,               &! fresh snow density perturbations
  alP,               &! albedo perturbations
  slP,               &! settling perturbations
  Sdird               ! Direct-beam shortwave radiation, per horizontal surface area (W/m2)
  
use GRID

use PARAMETERS

use PARAMMAPS

use SOILPARAMS 

use STATE_VARIABLES

use LANDUSE

implicit none
 
integer :: & 
  i,j,               &! Point counters
  k,                 &! Level counter
  iresults_count
  
integer :: &
  NNsmax,NNsoil,NNx,NNy
  
integer :: &
  NALBEDO,NCANMOD,NCONDCT,NDENSTY,NEXCHNG,NHYDROL,NSNFRAC,NRADSBG,NZOFFST,&
  NOSHDTN,NALRADT,NSNTRAN,NSNSLID,NSNOLAY,NCHECKS
  
real :: &
  DDs_min,DDs_surflay,zzT,zzU,zzRH

real :: &
  rtthresh

real :: &
  hcon_min            ! Thermal conductivity of soil minerals (W/m/K)
  
real, allocatable :: &
  fsat(:),           &! Initial moisture content of soil layers as fractions of saturation
  Tprof(:)            ! Initial soil layer temperatures (K)

character(len=200) :: nlst_file

character(len=20) :: CTILE
  
character(len=4), dimension(36) :: CLIST_DIAG_RESULTS

character(len=4), dimension(12) :: CLIST_STATE_RESULTS

real, dimension(:), allocatable :: &
  DDzsnow,DDzsoil

logical :: LHN_ON  ! activate the HN model

logical :: LZ0PERT,LWCPERT,LFSPERT,LALPERT,LSLPERT,LFOR_HN

logical :: lexist

!-1- !!!!!!!!!!!!!!!!!!!!  READ THE NAMELIST  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
namelist  /nam_grid/    NNx,NNy,NNsmax,NNsoil,DDs_min,DDs_surflay
namelist  /nam_layers/  DDzsnow,DDzsoil
namelist  /nam_driving/ zzT,zzU,zzRH
namelist  /nam_modconf/ NALBEDO,NCANMOD,NCONDCT,NDENSTY,NEXCHNG,NHYDROL,NSNFRAC,NRADSBG,NZOFFST,&
                        NOSHDTN,LHN_ON,LFOR_HN,NALRADT,NSNTRAN,NSNSLID,NSNOLAY,NCHECKS
namelist  /nam_modtile/ CTILE, rtthresh
namelist  /nam_modpert/ LZ0PERT,LWCPERT,LFSPERT,LALPERT,LSLPERT
namelist  /nam_results/ CLIST_DIAG_RESULTS, CLIST_STATE_RESULTS

! get namelist path from first command argument.
call getarg(1, nlst_file)
INQUIRE (FILE=nlst_file, EXIST=lexist)
if (lexist) then
  open(5000, file = nlst_file)
else
  print*, '  namelist file: ', trim(nlst_file), ' does not exist'
  call exit(1)
endif

! Initialization of variables
NNsmax = 3
NNsoil = 4 ! Attention, 5 soil layers in JIM ... 
NNx = 1
NNy = 1
DDs_min = 0.01 ! Minimum possible snow layer thickness (m)
DDs_surflay = 0.5 ! Maximum thickness of surface fine snow layering (m)
read(5000, nam_grid)
Nsmax = NNsmax
Nsoil = NNsoil
Nx = NNx
Ny = NNy
Ds_min = DDs_min
Ds_surflay = DDs_surflay

allocate(DDzsnow(Nsmax))
allocate(DDzsoil(Nsoil))
allocate(Dzsnow(Nsmax))
allocate(Dzsoil(Nsoil))
read(5000, nam_layers)
Dzsnow = DDzsnow
Dzsoil = DDzsoil

! Driving data
dt = 1
zzT = 10
zzU = 10
zzRH = 10
read(5000,nam_driving)
zT = zzT
zU = zzU
zRH = zzRH

! Model configuration
! -1 for mandatory NLST arguments.
NALBEDO = -1
NCANMOD = -1
NCONDCT = -1
NDENSTY = -1
NEXCHNG = -1
NHYDROL = -1
NSNFRAC = -1
NRADSBG = -1
NZOFFST = -1
NOSHDTN = -1
NALRADT = -1
NSNTRAN = -1
NSNSLID = -1
NSNOLAY = -1
NCHECKS = -1
LHN_ON = .FALSE.
LFOR_HN = .FALSE.
read(5000, nam_modconf)
ALBEDO = NALBEDO
CANMOD = NCANMOD
CONDCT = NCONDCT
DENSTY = NDENSTY
EXCHNG = NEXCHNG
HYDROL = NHYDROL
SNFRAC = NSNFRAC
RADSBG = NRADSBG
ZOFFST = NZOFFST
OSHDTN = NOSHDTN
ALRADT = NALRADT
SNTRAN = NSNTRAN
SNSLID = NSNSLID
SNOLAY = NSNOLAY
CHECKS = NCHECKS
HN_ON = LHN_ON
FOR_HN = LFOR_HN

if (ALBEDO==-1 .or. CANMOD==-1 .or. CONDCT==-1 .or. DENSTY==-1 .or. EXCHNG==-1 &
    .or. HYDROL==-1 .or. SNFRAC==-1 .or. RADSBG==-1 .or. ZOFFST==-1 .or. OSHDTN ==-1 .or. ALRADT ==-1 &
    .or. SNTRAN==-1 .or. SNSLID==-1 .or. SNOLAY==-1 .or. CHECKS ==-1) then
  print*, 'model configuration error:\n please specify all the fields of MODCONF in the namelist (&nam_modconf)'
  call exit(1)
endif

! Model perturbations
LZ0PERT = .FALSE.
LWCPERT = .FALSE.
LFSPERT = .FALSE.
LALPERT = .FALSE.
LSLPERT = .FALSE.
read(5000, nam_modpert)
Z0PERT = LZ0PERT
WCPERT = LWCPERT
FSPERT = LFSPERT
ALPERT = LALPERT
SLPERT = LSLPERT

! Modelled tile 
TILE = 'open'  ! set open as default if no tile has been defined 
tthresh = 0.1
read(5000, nam_modtile)
TILE = CTILE
tthresh = rtthresh

! List of output variables
CLIST_DIAG_RESULTS(:)  = '    ' ! BC 13-char length
CLIST_STATE_RESULTS(:) = '    ' ! BC 13-char length
read(5000, nam_results)
iresults_count = count(CLIST_DIAG_RESULTS /= '             ') ! BC 13-char length
LIST_DIAG_RESULTS(1:iresults_count) = CLIST_DIAG_RESULTS(1:iresults_count)
iresults_count = count(CLIST_STATE_RESULTS /= '             ') ! BC 13-char length
LIST_STATE_RESULTS(1:iresults_count) = CLIST_STATE_RESULTS(1:iresults_count)

! Outputs
Nave = 1 !24 Set to one, 24h averages created in OSHD Matlab wrapper

close(5000)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-2- !!!!!!!!!!!!!!!!!!!!    OPEN THE FILES   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Open Files for OSHD states input and output and results output
call OPEN_FILES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-3- !!!!!!!!!!!!!!!!!!!!  ALLOCATE VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(LW(Nx,Ny))
allocate(Ps(Nx,Ny))
allocate(Qa(Nx,Ny))
allocate(Rf(Nx,Ny))
allocate(Sf(Nx,Ny))
allocate(Sdir(Nx,Ny))
allocate(Sdif(Nx,Ny))
allocate(Ta(Nx,Ny))
allocate(Ua(Nx,Ny))
if (SNTRAN == 1) allocate(Udir(Nx,Ny))
allocate(RH(Nx,Ny))
allocate(Sf24h(Nx,Ny))
allocate(Tv(Nx,Ny))

if (Z0PERT) allocate(z0P(Nx,Ny))
if (WCPERT) allocate(wcP(Nx,Ny))
if (FSPERT) allocate(fsP(Nx,Ny))
if (ALPERT) allocate(alP(Nx,Ny))
if (SLPERT) allocate(slP(Nx,Ny))
if ((ALRADT == 1) .OR. (OSHDTN == 1)) allocate(Sdird(Nx,Ny))

! use Tv dummy in case of open simulations
if (CANMOD == 0) then
  Tv(:,:) = 1
endif

! Defaults for numerical solution parameters
Nitr = 4

! Defaults for canopy parameters
avg0 = 0.1
avgs = 0.4
cden = 0.004
cvai = 4.4
cveg = 20
Gcn1 = 0.5
Gcn2 = 0
gsnf = 0
kdif = 0.5
kveg = 1
rchd = 0.67
rchz = 0.2          
tcnc = 240
tcnm = 48

! Defaults for snow parameters
a_eta = 0.1
asmx = 0.86       ! unused if OSHDTN = 1
asmn = 0.6
b_eta = 0.023
bstb = 5
bthr = 2
c_eta = 250
eta0 = 3.7e7
eta1 = 7.62237e6
hfsn = 0.1
kfix = 0.24
rho0 = 300
rhob = 6
rhoc = 26
rhof = 109
rhos_min = 50.0
rhos_max = 750.0
rcld = 300
rgr0 = 5e-5
rmlt = 500
snda = 2.8e-6
Talb = -2
tcld = 1000
tmlt = 100
trho = 200
Wirr = 0.03
z0sn = 0.002
Sfmin = 10

! some defaults different for forest tile - commented-out values based on FS-EBP runs, revisit during tuning
if (TILE == 'forest') then
  ! asmx = 0.88
  hfsn = 0.3    
  z0sn = 0.01   
endif 

! Defaults for ground surface parameters
bstb = 5
gsat = 0.01

! Defaults for additional parameters required for forest snow process parametrization
adfs = 3
adfl = 2
fsar = 0.1
psf  = 1
psr  = 0.1
wcan = 2.5
zsub = 2
zgf = 1
zgr = 0
khcf = 3

if (DENSTY == 0) then
  rhof = rho0
endif

! Surface data from defaults or namelists
! Surface properties 
allocate(alb0(Nx,Ny))
allocate(fcly(Nx,Ny))
allocate(fsnd(Nx,Ny))
allocate(z0sf(Nx,Ny))
if (SNTRAN == 1) allocate(vegsnowd_xy(Nx,Ny))
if (TILE == 'glacier') then
  alb0(:,:) = 0.3
  z0sf(:,:) = 0.04
else
  alb0(:,:) = 0.2
  z0sf(:,:) = 0.2
endif
fcly(:,:) = 0.3
fsnd(:,:) = 0.6
if (TILE == 'forest') z0sf(:,:) = 0.2
if (SNTRAN == 1) vegsnowd_xy(:,:) = 0.1

! Canopy parameters
allocate(canh(Nx,Ny))
allocate(fsky(Nx,Ny))
allocate(fveg(Nx,Ny)) 
allocate(fves(Nx,Ny))
allocate(hcan(Nx,Ny))
allocate(lai(Nx,Ny))
allocate(pmultf(Nx,Ny))
allocate(scap(Nx,Ny))
allocate(trcn(Nx,Ny))
allocate(VAI(Nx,Ny))
allocate(vfhp(Nx,Ny))
canh(:,:)        = undef
fsky(:,:)        = undef
fveg(:,:)        = undef
fves(:,:)        = undef
hcan(:,:)        = undef
pmultf(:,:)      = undef
scap(:,:)        = undef
trcn(:,:)        = undef
VAI(:,:)         = undef

! Terrain properties
allocate(fsky_terr(Nx,Ny))
allocate(lat(Nx,Ny))
allocate(lon(Nx,Ny))
allocate(dem(Nx,Ny))
allocate(tilefrac(Nx,Ny))
allocate(glacierfrac(Nx,Ny))
fsky_terr(:,:)   = undef
lat(:,:)         = undef
lon(:,:)         = undef
dem(:,:)         = undef
tilefrac(:,:)    = undef
glacierfrac(:,:) = undef

if (SNFRAC == 0) then
  allocate(slopemu(Nx,Ny))
  allocate(xi(Nx,Ny))
  slopemu(:,:)   = undef
  xi(:,:)        = undef
endif
if (SNFRAC == 0 .or. SNTRAN == 1) then
  allocate(Ld(Nx,Ny))
  Ld(:,:)        = undef
endif

if (SNSLID == 1) then
  allocate(slope(Nx,Ny))
  allocate(Shd(Nx,Ny))
  slope(:,:)     = undef
  Shd(:,:)       = undef
endif

! Derived soil parameters
allocate(b(Nx,Ny))
allocate(hcap_soil(Nx,Ny))
allocate(hcon_soil(Nx,Ny))
allocate(sathh(Nx,Ny))
allocate(Vsat(Nx,Ny))
allocate(Vcrit(Nx,Ny))
do j = 1, Ny
  do i = 1, Nx
    if (fcly(i,j) + fsnd(i,j) > 1) then
      fcly(i,j) = 1 - fsnd(i,j)
    endif
    b(i,j) = 3.1 + 15.7*fcly(i,j) - 0.3*fsnd(i,j)
    hcap_soil(i,j) = (2.128*fcly(i,j) + 2.385*fsnd(i,j))*1e6 / (fcly(i,j) + fsnd(i,j))
    sathh(i,j) = 10**(0.17 - 0.63*fcly(i,j) - 1.58*fsnd(i,j))
    Vsat(i,j) = 0.505 - 0.037*fcly(i,j) - 0.142*fsnd(i,j)
    Vcrit(i,j) = Vsat(i,j)*(sathh(i,j)/3.364)**(1/b(i,j))
    hcon_min = (hcon_clay**fcly(i,j)) * (hcon_sand**(1 - fcly(i,j)))
    hcon_soil(i,j) = (hcon_air**Vsat(i,j)) * (hcon_min**(1 - Vsat(i,j)))
  end do
end do

! Convert time scales from hours to seconds
dt = 3600*dt
tcnc = 3600*tcnc
tcnm = 3600*tcnm
tcld = 3600*tcld
tmlt = 3600*tmlt
trho = 3600*trho

! Allocate state variables
allocate(albs(Nx,Ny))
allocate(Ds(Nsmax,Nx,Ny))
allocate(Nsnow(Nx,Ny))
allocate(Qcan(Nx,Ny))
allocate(rgrn(Nsmax,Nx,Ny))
allocate(histowet(Nsmax,Nx,Ny)) ! LQ: histowet allocated in any case to simplify the layering routine
allocate(Sice(Nsmax,Nx,Ny))
allocate(Sliq(Nsmax,Nx,Ny))
allocate(Sveg(Nx,Ny))
allocate(Tcan(Nx,Ny))
allocate(theta(Nsoil,Nx,Ny))
allocate(Tsnow(Nsmax,Nx,Ny))
allocate(Tsoil(Nsoil,Nx,Ny))
allocate(Tsrf(Nx,Ny))
allocate(fsnow(Nx,Ny))
allocate(Tveg(Nx,Ny))

! Default initialization of state variables 
albs(:,:)               = undef
Ds(:,:,:)               = undef
fsnow(:,:)              = undef
Nsnow(:,:)              = iundef
Qcan(:,:)               = undef
rgrn(:,:,:)             = undef !*GM watch out: rgrn currently not tracked
histowet(:,:,:)         = undef ! LQ: histowet only tracked if SNTRAN
Sice(:,:,:)             = undef
Sliq(:,:,:)             = undef
Sveg(:,:)               = undef
Tcan(:,:)               = undef
Tsnow(:,:,:)            = undef
Tsoil(:,:,:)            = undef
Tveg(:,:)               = undef

! Allocation and initialization of optional state variables
if (SNFRAC == 0 .or. SNFRAC == 2) then
  allocate(snowdepthmax(Nx,Ny))
  snowdepthmax(:,:) = undef
endif

if (SNFRAC == 0) then
  allocate(snowdepthmin(Nx,Ny))
  allocate(snowdepthhist(14,Nx,Ny))
  allocate(swemin(Nx,Ny))
  allocate(swemax(Nx,Ny))
  allocate(swehist(14,Nx,Ny))
  snowdepthmin(:,:)     = undef
  snowdepthhist(:,:,:)  = undef
  swemin(:,:)           = undef
  swemax(:,:)           = undef
  swehist(:,:,:)        = undef
endif

if (SNTRAN == 1) then
  allocate(dSWE_tot_subl(Nx,Ny))
  allocate(dSWE_tot_salt(Nx,Ny))
  allocate(dSWE_tot_susp(Nx,Ny))
  dSWE_tot_subl(:,:)    = undef
  dSWE_tot_salt(:,:)    = undef
  dSWE_tot_susp(:,:)    = undef
endif

if (SNSLID == 1) then
  allocate(dSWE_tot_slide(Nx,Ny))
  allocate(index_sorted_dem(Nx*Ny,2))
  dSWE_tot_slide(:,:)   = undef
  index_sorted_dem(:,:) = iundef
endif

! Initial soil profiles from namelist
allocate(fsat(Nsoil))
allocate(Tprof(Nsoil))
fsat(:)  = 0.5
Tprof(:) = 285
do k = 1, Nsoil
  theta(k,:,:) = fsat(k)*Vsat(:,:)
  Tsoil(k,:,:) = Tprof(k)
end do
Tsrf(:,:) = Tsoil(1,:,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-4- !!!!!!!!!!!!!!!!!!!! READ DRIVING/STATES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! states relevant to both open and forest simulation
read(1101) ((albs(i,j),i=1,Nx),j=1,Ny)
read(1102) (((Ds(k,i,j),k=1,Nsmax),i=1,Nx),j=1,Ny)
read(1103) ((fsnow(i,j),i=1,Nx),j=1,Ny)
read(1104) ((Nsnow(i,j),i=1,Nx),j=1,Ny)
read(1106) (((Sice(k,i,j),k=1,Nsmax),i=1,Nx),j=1,Ny)
read(1107) (((Sliq(k,i,j),k=1,Nsmax),i=1,Nx),j=1,Ny)
read(1116) ((Tsrf(i,j),i=1,Nx),j=1,Ny)
read(1119) (((Tsnow(k,i,j),k=1,Nsmax),i=1,Nx),j=1,Ny)
read(1120) (((Tsoil(k,i,j),k=1,Nsoil),i=1,Nx),j=1,Ny)
read(1123) ((fsky_terr(i,j),i=1,Nx),j=1,Ny)
read(1127) ((lat(i,j),i=1,Nx),j=1,Ny)
read(1128) ((lon(i,j),i=1,Nx),j=1,Ny)
read(1129) ((dem(i,j),i=1,Nx),j=1,Ny)

! Cap glacier temperatures to 0°C 
if (TILE == 'glacier') then
  do j = 1, Ny
  do i = 1, Nx
    Tsrf(i,j) = min(Tsrf(i,j),Tm)
    do k = 1, Nsoil
      Tsoil(k,i,j) = min(Tsoil(k,i,j),Tm)
    end do
  end do
  end do
endif

! model tile fractions 
if (TILE == 'open') then 
  tilefrac = dem/dem   ! temporary fix to get ones within our entire domain, assuming we always want to run an open tile. may have to be revisited
  if (SNTRAN == 1 .OR. SNSLID == 1) then
    read(1148) ((glacierfrac(i,j),i=1,Nx),j=1,Ny)
  endif
else 
  read(1139) ((tilefrac(i,j),i=1,Nx),j=1,Ny)
endif 

if (TILE == 'open' .AND. (SNTRAN == 1 .OR. SNSLID == 1)) then ! Cap glacier temperatures only in relevant pixels, and alter snow free albedo and roughness lengths
  do j = 1, Ny
  do i = 1, Nx
    if (glacierfrac(i,j) > epsilon(glacierfrac(i,j))) then
      Tsrf(i,j) = min(Tsrf(i,j),Tm)
      alb0(i,j) = 0.3 
      z0sf(i,j) = 0.04
      do k = 1, Nsoil
        Tsoil(k,i,j) = min(Tsoil(k,i,j),Tm)
      end do
    endif
  end do
  end do
endif 

if (SNFRAC == 0 .or. SNFRAC == 2) then
  read(1110) ((snowdepthmax(i,j),i=1,Nx),j=1,Ny)
endif

if (SNFRAC == 0) then
  ! states specific to open runs
  read(1109) ((snowdepthmin(i,j),i=1,Nx),j=1,Ny)
  read(1111) (((snowdepthhist(k,i,j),k=1,14),i=1,Nx),j=1,Ny)
  read(1113) ((swemin(i,j),i=1,Nx),j=1,Ny)
  read(1114) ((swemax(i,j),i=1,Nx),j=1,Ny)
  read(1115) (((swehist(k,i,j),k=1,14),i=1,Nx),j=1,Ny)
  read(1124) ((slopemu(i,j),i=1,Nx),j=1,Ny)
  read(1125) ((xi(i,j),i=1,Nx),j=1,Ny)
endif

if (SNFRAC == 0 .or. SNTRAN == 1) then
  read(1126) ((Ld(i,j),i=1,Nx),j=1,Ny)
endif

if (TILE /= 'forest') then
  ! canopy properties (no canopy)
  VAI(:,:)  = 0
  hcan(:,:) = 0
  fsky(:,:) = 1
  trcn(:,:) = exp(-kdif*VAI(:,:))
  fveg(:,:) = 1 - exp(-kveg*VAI(:,:))
  fves(:,:) = 1 - exp(-kveg*VAI(:,:))
else ! TILE == 'forest'
  ! lus fields specific to forest runs
  read(1130) ((Qcan(i,j),i=1,Nx),j=1,Ny)
  read(1131) ((Sveg(i,j),i=1,Nx),j=1,Ny)
  read(1132) ((Tcan(i,j),i=1,Nx),j=1,Ny)
  read(1133) ((Tveg(i,j),i=1,Nx),j=1,Ny)
  read(1134) ((fveg(i,j),i=1,Nx),j=1,Ny)
  read(1135) ((hcan(i,j),i=1,Nx),j=1,Ny)
  read(1136) ((lai(i,j),i=1,Nx),j=1,Ny)
  read(1137) ((vfhp(i,j),i=1,Nx),j=1,Ny)
  read(1138) ((fves(i,j),i=1,Nx),j=1,Ny)
  read(1140) ((pmultf(i,j),i=1,Nx),j=1,Ny)

  ! derived canopy properties 
  VAI(:,:) = lai(:,:) 
  trcn(:,:) = 1-0.9*fveg(:,:)  
  do j = 1, Ny
    do i = 1, Nx
      fsky(i,j) = vfhp(i,j)/trcn(i,j)
      if ( fsky(i,j) > 1 ) trcn(i,j) = vfhp(i,j)
      if ( fsky(i,j) > 1 ) fsky(i,j) = 1
    end do
  end do 
endif

! derived canopy parameters
canh(:,:) = 12500*VAI(:,:)
scap(:,:) = cvai*VAI(:,:)

if (SNTRAN == 1) then
  ! states specific to SNOWTRAN3D
  read(1121) (((histowet(k,i,j),k=1,Nsmax),i=1,Nx),j=1,Ny)
  read(1141) ((dSWE_tot_subl(i,j),i=1,Nx),j=1,Ny)
  read(1142) ((dSWE_tot_salt(i,j),i=1,Nx),j=1,Ny)
  read(1143) ((dSWE_tot_susp(i,j),i=1,Nx),j=1,Ny)
endif

if (SNSLID == 1) then
  ! states specific to SnowSlide
  read(1144) ((slope(i,j),i=1,Nx),j=1,Ny)
  read(1145) ((Shd(i,j),i=1,Nx),j=1,Ny)
  read(1146) ((dSWE_tot_slide(i,j),i=1,Nx),j=1,Ny)
  read(1147) ((index_sorted_dem(i,k),i=1,Nx*Ny),k=1,2)
endif


! Tuned snow surface properties

allocate(adc(Nx,Ny))
allocate(afs(Nx,Ny))
allocate(z0_snow(Nx,Ny))

if (OSHDTN == 0) then

  adm = 100
  adc(:,:) = 1000
  afs(:,:) = asmx
  if (TILE == 'glacier' .or. ((SNTRAN == 1 .or. SNSLID == 1) .and. glacierfrac(i,j) > epsilon(glacierfrac(i,j)))) then
    z0_snow(:,:) = 0.0009
  else
    z0_snow(:,:) = z0sn
  endif

else ! OSHDTN == 1

  adm = 130

  do j = 1, Ny
    do i = 1, Nx
    
      ! Elevation-dependent tuning of cold snow albedo decay time
      if (dem(i,j) >= 2300) then
        adc(i,j)  = 6000
      elseif (dem(i,j) <= 1500) then
        adc(i,j) = 3000
      else
        adc(i,j) = 6000 + (2300 - dem(i,j)) / (2300 - 1500) * (3000 - 6000)
      endif
      
      ! Fresh snow albedo is now constant (previously elevation-dependent)
      afs(i,j) = asmx
      
      ! Elevation-dependent tuning of snow roughness length
      if (TILE == 'glacier' .or. ((SNTRAN == 1 .or. SNSLID == 1) .and. glacierfrac(i,j) > epsilon(glacierfrac(i,j)))) then
        z0_snow(i,j) = 0.0009
      elseif (TILE == 'forest') then
        z0_snow(i,j) = z0sn
      else
        if (dem(i,j) >= 2300) then
          z0_snow(i,j) = 0.01
        elseif (dem(i,j) >= 1500) then
          z0_snow(i,j) = 0.2 + (dem(i,j) - 1500) / (2300 - 1500) * (0.01 - 0.2)
        else
          z0_snow(i,j) = 0.2
        endif
      endif
      
    end do
  end do
  
endif

end subroutine SETUP

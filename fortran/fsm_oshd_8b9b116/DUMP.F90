!-----------------------------------------------------------------------
! Write out state variables at end of run
!-----------------------------------------------------------------------
subroutine DUMP

use MODCONF, only: CANMOD, SNFRAC, SNTRAN, SNSLID

use MODE_WRITE, only: WRITE_2D

use STATE_VARIABLES, only : &
  albs,                    &! Snow albedo
  Ds,                      &! Snow layer thicknesses (m)
  histowet,                &! Historical variable for past wetting of a layer (0-1)
  fsnow,                   &! Snow cover fraction
  Nsnow,                   &! Number of snow layers 
  Qcan,                    &! Canopy air space humidity
  Sice,                    &! Ice content of snow layers (kg/m^2)
  Sliq,                    &! Liquid content of snow layers (kg/m^2)
  snowdepthmin,            &! min Snow at time step of swemin (m)
  snowdepthmax,            &! max Snow at time step of swemax (m)
  snowdepthhist,           &! history of Snow depth during last 14 days (m)  
  swemin,                  &! min swe during season (mm)
  swemax,                  &! max swe during season (mm)  
  swehist,                 &! history of Snow depth during last 14 days (kg/m^2)
  Sveg,                    &! Canopy snow mass (kg/m^2)
  Tcan,                    &! Canopy air space temperature (K)
  Tsnow,                   &! Snow layer temperatures (K)
  Tsoil,                   &! Soil layer temperatures (K)
  Tsrf,                    &! Surface skin temperature (K)
  Tveg,                    &! Vegetation temperature (K)
  dSWE_tot_subl,           &! Cumulated SWE change due to sublimation (kg/m^2)
  dSWE_tot_salt,           &! Cumulated SWE change due to saltation (kg/m^2)
  dSWE_tot_susp,           &! Cumulated SWE change due to suspension (kg/m^2)
  dSWE_tot_slide,          &! Cumulated SWE change due to snow slides (kg/m^2)
  index_sorted_dem          ! Location (i,j) of sorted grid points

use GRID, only: &
  Nsmax,                   &! Maximum number of snow layers
  Nsoil,                   &! Number of soil layers
  Nx,Ny                     ! Grid dimensions

implicit none

integer :: i, j, k
  
! Write into state files.
write(1201) ((albs(i,j),i=1,Nx),j=1,Ny)
write(1202) (((Ds(k,i,j),k=1,Nsmax),i=1,Nx),j=1,Ny)
write(1203) ((fsnow(i,j),i=1,Nx),j=1,Ny)
write(1204) ((Nsnow(i,j),i=1,Nx),j=1,Ny)
write(1206) (((Sice(k,i,j),k=1,Nsmax),i=1,Nx),j=1,Ny)
write(1207) (((Sliq(k,i,j),k=1,Nsmax),i=1,Nx),j=1,Ny)
write(1216) ((Tsrf(i,j),i=1,Nx),j=1,Ny)
write(1219) (((Tsnow(k,i,j),k=1,Nsmax),i=1,Nx),j=1,Ny)
write(1220) (((Tsoil(k,i,j),k=1,Nsoil),i=1,Nx),j=1,Ny)
if (SNFRAC == 0 .or. SNFRAC == 2) then
  write(1210) ((snowdepthmax(i,j),i=1,Nx),j=1,Ny)
endif
if (SNFRAC == 0) then
  write(1209) ((snowdepthmin(i,j),i=1,Nx),j=1,Ny)
  write(1211) (((snowdepthhist(k,i,j),k=1,14),i=1,Nx),j=1,Ny)
  write(1213) ((swemin(i,j),i=1,Nx),j=1,Ny)
  write(1214) ((swemax(i,j),i=1,Nx),j=1,Ny)
  write(1215) (((swehist(k,i,j),k=1,14),i=1,Nx),j=1,Ny)
endif
if (CANMOD == 1) then
  write(1223) ((Qcan(i,j),i=1,Nx),j=1,Ny)
  write(1224) ((Sveg(i,j),i=1,Nx),j=1,Ny)
  write(1225) ((Tcan(i,j),i=1,Nx),j=1,Ny)
  write(1226) ((Tveg(i,j),i=1,Nx),j=1,Ny)
endif
if (SNTRAN == 1) then
  ! states specific to SNOWTRAN3D
  write(1221) (((histowet(k,i,j),k=1,Nsmax),i=1,Nx),j=1,Ny)
  write(1250) ((dSWE_tot_subl(i,j),i=1,Nx),j=1,Ny)
  write(1251) ((dSWE_tot_salt(i,j),i=1,Nx),j=1,Ny)
  write(1252) ((dSWE_tot_susp(i,j),i=1,Nx),j=1,Ny)
endif
if (SNSLID == 1) then
  ! states specific to SnowSlide
  write(1253) ((dSWE_tot_slide(i,j),i=1,Nx),j=1,Ny)
  write(1254) ((index_sorted_dem(i,k),i=1,Nx*Ny),k=1,2)
endif

end subroutine DUMP

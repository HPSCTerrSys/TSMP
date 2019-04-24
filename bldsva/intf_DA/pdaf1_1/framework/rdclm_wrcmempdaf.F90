module rdclm_wrcmempdaf 

  use netcdf
  use clm4cmem, only: CLM_DATA, SHOWINFO
  use parkind1, only: JPIM, JPRM
  USE YOMCMEMNETCDF, ONLY : NTIMES_SM,NLATS_SM, NLONS_SM, NLVLS_SM,NTIMES, xlats,xlons,xlvls,xtimes
contains

  ! ------------------------------------------------------------------
  ! SUBROUTINE TO READ CLM INPUT FILE from PDAF memory AND PROCESS THE DATA FOR CMEM
  ! ------------------------------------------------------------------
  ! INPUTS:
  ! - LS   : state vector for assimilation
  ! - SAT  : satellites information
  ! - i_idx : gridcell index
  ! OUTPUTS: 
  ! - LS        : the structure contains all inputs
  !-------------------------------------------------------------------
  
  subroutine read_CLM_pdaf(LS,SAT,Nvars)

!#if model == tag_model_clm
   USE PARKIND1,      ONLY : JPIM
   USE clmtype      , only : clm3
   USE clm_varpar   , only : nlevsoi, nlevlak,nlevsno
   use shr_kind_mod,  only: r8 => shr_kind_r8
   use clm4cmem,      only: SATELLITE, local_incidence, gradientm
   use domainMod,     only: ldomain
   use decompMod,     only: get_proc_bounds
   ! use surfrdMod, only: pctpft
   use mod_tsmp,      only: idx_map_subvec2state_fortran,pf_statevecsize
   use mod_assimilation,  ONLY: obs_index_p
   implicit none
   type(SATELLITE), intent(inout) :: SAT
   type(CLM_DATA), intent(out) :: LS
   integer(kind=JPIM), intent(in) :: Nvars(3)


   integer :: begp, endp   ! per-proc beginning and ending pft indices
   integer :: begc, endc   ! per-proc beginning and ending column indices
   integer :: begl, endl   ! per-proc beginning and ending landunit indices
   integer :: begg, endg   ! per-proc gridcell ending gridcell indices
   integer :: clm_begg,clm_endg 
   integer :: nlev, idxtime=1
   
   integer :: NINC, NLATS, NLONS, NTIME, NLEVPFT
   real(r8), pointer :: latdeg(:)       !latitude (degrees)
   real(r8), pointer :: londeg(:)       !longitude (degrees)
   real(r8), pointer :: dz(:,:)         !layer thickness (m)  (-nlevsno+1:nlevgrnd) 
   real(r8), pointer :: z(:,:)          !layer depth (m) (-nlevsno+1:nlevgrnd) 
   !real(r8), pointer :: dz_lake(:,:)    !layer thickness (m)  (-nlevsno+1:nlevgrnd)
   !real(r8), pointer :: z_lake(:,:)     ! layer depth for lake (m)
   real(r8), pointer :: swc(:,:)
   real(r8), pointer :: psand(:,:)
   real(r8), pointer :: pclay(:,:)
   real(r8), pointer :: snowdp(:)
   real(r8), pointer :: snowliq(:)
   real(r8), pointer :: t_soisno(:,:)
   real(r8), pointer :: t_grnd(:)
   real(r8), pointer :: tlai(:)        
   real(r8), pointer :: t_ref2m(:)
   integer , pointer :: itype(:)
   real(r8), pointer :: topo(:)
   ! real(r8), pointer :: longxy(:)
   ! real(r8), pointer :: latixy(:)
   real(r8), pointer :: time(:)
   real(r8), pointer :: fcov(:) 
   integer  :: i,j,cc=1,offset=0
   real(r8) :: swcave
   real(r8) :: swctmp(10)

   !real, dimension(:), allocatable :: lats, lons, time, levlak, levsoi
   !real, dimension(:,:), allocatable :: Z, WATER, LSM, CVL, CVH, TVL, TVH, CLAY, SAND   ! 2D var
   !real, dimension(:,:), allocatable :: latixy, longxy, slope, aspect
   !real, dimension(:,:,:), allocatable :: TMP, LAI, theta_inc  ! 3D var
   !real, dimension(:,:,:,:), allocatable :: TMP4D, LAIL     ! 4D var
   !real, dimension(:,:,:,:), allocatable :: SD, RSN, STL  ! 4D var
   !real, allocatable :: TSKIN(:,:,:,:),SWVL(:,:,:,:), TAIR(:,:,:,:)  ! 4D var
   !level indexes e.g. idxlev(nlev) = (/1, 2, 3/) 1st, 2nd, 3rd levels
   real, dimension(:), allocatable :: levsoi,levlak
   integer, dimension(:), allocatable :: idxlev
   
   latdeg   => clm3%g%latdeg 
   londeg   => clm3%g%londeg 
   dz       => clm3%g%l%c%cps%dz
   !dz_lake  => clm3%g%l%c%cps%dz_lake
   z        => clm3%g%l%c%cps%z
   !z_lake   => clm3%g%l%c%cps%z_lake
   swc      => clm3%g%l%c%cws%h2osoi_vol
   psand    => clm3%g%l%c%cps%psand
   pclay    => clm3%g%l%c%cps%pclay
   snowdp   => clm3%g%l%c%cps%snowdp
   snowliq  => clm3%g%l%c%cws%snowliq
   t_soisno => clm3%g%l%c%ces%t_soisno
   t_grnd   => clm3%g%l%c%ces%t_grnd
   tlai     => clm3%g%l%c%p%pps%tlai
   t_ref2m  => clm3%g%l%c%p%pes%t_ref2m
   itype    => clm3%g%l%c%itype 
   time     => clm3%g%l%c%p%pepv%days_active
   fcov     => clm3%g%l%c%cws%fcov         
   topo     => ldomain%topo
   ! longxy   => ldomain%longxy
   ! latixy   => ldomain%latixy

   
   ! call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
   clm_begg     = 1
   clm_endg     = Nvars(1)*Nvars(2)  
   NLATS        = Nvars(2)
   NLONS        = Nvars(1)
   NINC         = 1 
   NTIME        = 1
   NLEVPFT      = 1 
   print*,  "clm_begg/clm_endg is ",clm_begg,clm_endg

   ! Get the varids of the latitude and longitude coordinate variables.
   ! Read the latitude and longitude data.
   allocate(LS%lats(NLATS))
   allocate(LS%lons(NLONS))
   do j=1,NLATS
      LS%lats(j) = latdeg((j-1)*NLONS+1)
   end do 
   do j=1,NLONS
      LS%lons(j) = londeg(j)
   end do 
   print*, 'input LS LAT/LON begins at ',LS%lats(1),LS%lons(1)
   print*, 'input LS LAT/LON ends at ',LS%lats(NLATS),LS%lons(NLONS)

   ! Setting default values for input variables in case not passed in call  
   nlev         = nlevsoi
   print*, 'input LS NLEV is ',nlev
   allocate(idxlev(nlev))
   idxlev = (/(I,I=1,nlev)/)
   if(.not.allocated(LS%levsoi)) allocate(LS%levsoi(nlev))
   LS%levsoi = z(1,idxlev)
   !if(.not.allocated(LS%levlak)) allocate(LS%levlak(nlev))
   !LS%levlak = dz_lake(1,idxlev)
   !print*, 'input LS levlak is ',LS%levlak
   LS%DATESTRING = '18/04/20 10:28:22'
   if(.not.allocated(LS%time)) allocate(LS%time(NTIME))
   LS%time = time(idxtime)
   NINC = size(SAT%theta)
   if(allocated(SAT%time)) deallocate(SAT%time)
   allocate(SAT%time(NTIME))
   SAT%time = time(idxtime)
   print*, 'SAT time is ',SAT%time 
   
   ! 0 -> Reading TOPO:
   allocate(LS%Z(NLONS,NLATS))
   ! geopotential height from m(CLM) to km(CMEM)
   do i = 1,NLONS
      do j = 1,NLATS
         LS%Z(i,j) = topo((j-1)*NLONS+i)/1000
      enddo
   enddo
   print*, 'geopotential height(km) is ',LS%Z(1,1),LS%Z(NLONS,NLATS)

   ! 0.1 -> reading the longitude and latitude grid-boxes:
   allocate(LS%latixy(NLONS,NLATS), LS%longxy(NLONS,NLATS))
   LS%latixy = reshape(latdeg, (/NLONS,NLATS/))
   LS%longxy = reshape(londeg, (/NLONS,NLATS/))
   print*, 'reading the longitude and latitude is done'
   print*, LS%latixy(1,1),LS%latixy(NLONS,NLATS)   
   print*, LS%longxy(1,1),LS%longxy(NLONS,NLATS)
   ! 0.2 -> calculating slope and aspect from TOPO Z:
   allocate(LS%slope(NLONS,NLATS), LS%aspect(NLONS,NLATS))
   allocate(LS%theta_inc(NLONS,NLATS,NINC))
   call gradientm(LS%Z, LS%longxy, LS%latixy, LS%slope,&
        LS%aspect, LS%resol_km)

   ! 0.3 -> calculating pixel-wise incidence angle theta_inc:
   LS%theta = SAT%theta
   LS%theta_inc = local_incidence(LS%slope, LS%aspect, SAT)
   print*, 'theta is ',LS%theta
   ! 1. Making up soil moisture [%/%]        
   allocate(LS%SWVL(NLONS,NLATS,NTIME,nlev))
   LS%SWVL = reshape(swc, (/NLONS,NLATS,NTIME,nlev/))
   print*, 'soil moisture is ',LS%SWVL(1,1,1,1),LS%SWVL(NLONS,NLATS,NTIME,nlev)
   ! 2. Making up soil temperature [K]
   allocate(LS%STL(NLONS,NLATS,NTIME,nlev))  
   !LS%STL = reshape(t_soisno(:,nlevsno+1:nlevsno+nlev),(/NLONS,NLATS,NTIME,nlev/))
   LS%STL = reshape(t_soisno(:,1:nlev),(/NLONS,NLATS,NTIME,nlev/))
   print*, 'soil temperature is',LS%STL(1,1,NTIME,nlev),LS%STL(NLONS,NLATS,NTIME,nlev)
   ! 3. Making up snow depth
   allocate(LS%SD(NLONS,NLATS,1,NTIME))
   LS%SD = reshape(snowdp, (/NLONS,NLATS,1,NTIME/)) 
   print*, 'snow depth is ',LS%SD(1,1,1,NTIME),LS%SD(NLONS,NLATS,1,NTIME)
   ! 4. Making up snow density [kg/m^3]
   !  water equivalent of snow/snow depth [m]
   allocate(LS%RSN(NLONS,NLATS,1,NTIME))
   do i = 1,NLONS
      do j = 1,NLATS 
         if (snowdp((j-1)*NLONS+i) .eq. 0) then
            LS%RSN(i,j,1,NTIME) = 0
         else
            LS%RSN(i,j,1,NTIME) = snowliq((j-1)*NLONS+i)/snowdp((j-1)*NLONS+i)
         end if
      end do
   end do
   print*, 'snow RSN is ',LS%RSN(1,1,1,NTIME),LS%RSN(NLONS,NLATS,1,NTIME)
   ! 5. Making up skin temperature
   allocate(LS%TSKIN(NLONS,NLATS,1,NTIME))
   LS%TSKIN = reshape(t_grnd, (/NLONS,NLATS,1,NTIME/))
   print*, 'skin temperature is',LS%TSKIN(1,1,1,NTIME),LS%TSKIN(NLONS,NLATS,1,NTIME) 
   ! 6 texture
   allocate(LS%CLAY(NLONS,NLATS))
   allocate(LS%SAND(NLONS,NLATS))
   do i  = 1,NLONS
      do j = 1,NLATS
         LS%SAND(i,j) = psand((i-1)*NLONS+j,1)
         LS%CLAY(i,j) = pclay((i-1)*NLONS+j,1)
      enddo
   enddo
   print*, 'sand fraction is ',LS%SAND(1,1),LS%SAND(NLONS,NLATS)
   print*, 'clay fraction is ',LS%CLAY(1,1),LS%CLAY(NLONS,NLATS)
   ! Vegetation
   ! Vegetation types    | CLM VR0  | ECOCLIMAP | TESSEL :
   ! ---------------------------------------------------------------
   ! bare_soil           | 0        |           |
   ! needle_leaf_forest  | 1        |  2        | 3      (High Veg)
   ! broad_leaf_forest   | 7        |  1        | 4      (High Veg)
   ! grassland           | 14       |  5        | 2      (Low Veg)
   ! crops/agricul_land  | 15       |  6        | 1      (Low Veg)
   ! ----------------------------------------------e-----------------
   ! 7. vegetation
   ! 7.1 ! 7.1 -> high/low vegetation fraction. TerrsysMP is high spatial
   ! resolution, so one pft for each cell. Set fracation to 1.  
  
   allocate(LS%CVL(NLONS,NLATS))   ! TESSEL fraction of low veg
   allocate(LS%CVH(NLONS,NLATS))   ! TESSEL fraction of high veg
   
   do i  = 1,NLONS
      do j = 1,NLATS
            
         if(itype((j-1)*NLONS+i).eq.1 .or. itype((j-1)*NLONS+i).eq.7) then
            LS%CVH(i,j) = 1
            LS%CVL(i,j) = 0
         elseif (itype((j-1)*NLONS+i).eq.14 .or. itype((j-1)*NLONS+i).eq.15) then
            LS%CVH(i,j) = 0
            LS%CVL(i,j) = 1
         else
            LS%CVH(i,j) = 0
            LS%CVL(i,j) = 0
         endif
      end do
   end do
    
   ! 7.2 high/low vegetation type
   allocate(LS%TVL(NLONS,NLATS))   ! type of low veg
   allocate(LS%TVH(NLONS,NLATS))   ! type of high veg 
   do i  = 1,NLONS       
      do j = 1,NLATS
       ! ECOCLIM: 2=Needleleaf, 1=Deciduous
       ! high vegeation first       
       ! index CLM 1 and 7
       ! ECOCLIM: 5=short_grass, 6=Crops
       ! low vegeation next
       ! index CLM 14 and 15
         if(itype((j-1)*NLONS+i).eq.2) then
            LS%TVH(i,j)   = 2
            LS%TVL(i,j)   = 0
         elseif (itype((j-1)*NLONS+i).eq.8) then
            LS%TVH(i,j)   = 1
            LS%TVL(i,j)   = 0
         elseif (itype((j-1)*NLONS+i).eq.15) then
            LS%TVH(i,j)   = 0
            LS%TVL(i,j)   = 5
         elseif (itype((j-1)*NLONS+i).eq.16) then
            LS%TVH(i,j)   = 0
            LS%TVL(i,j)   = 6
         else
            LS%TVH(i,j)   = 0
            LS%TVL(i,j)   = 0
         endif
      end do
   end do
   ! 7.3 water fraction [-] / land franction [-]
   if(allocated(LS%WATER)) deallocate(LS%WATER)
   allocate(LS%WATER(NLONS,NLATS))
   if(allocated(LS%LSM)) deallocate(LS%LSM)
   allocate(LS%LSM(NLONS,NLATS)) 
   do i = 1,NLONS 
      do j = 1, NLATS
      if(fcov((j-1)*NLONS+i).eq.1) then
            LS%WATER(i,j) = 1
            LS%LSM(i,j)   = 0
         else
            LS%WATER(i,j) = 0
            LS%LSM(i,j)   = 1
         endif
      end do
   end do

   ! 7.4 LAI of low vegetation
   if(allocated(LS%LAIL)) deallocate(LS%LAIL)
   allocate(LS%LAIL(NLONS,NLATS,1,NTIME))
   do i = 1,NLONS
      do j = 1, NLATS
         LS%LAIL(i,j,1,NTIME) = tlai((j-1)*NLONS+i)
      enddo
   enddo

   print*, 'LAI of low vegeation is ',LS%LAIL(1,1,1,NTIME),LS%LAIL(NLONS,NLATS,1,NTIME)
   ! 8. 2m Air Temperature [K]
   if (allocated(LS%TAIR)) deallocate(LS%TAIR)
   allocate(LS%TAIR(NLONS,NLATS,1,NTIME))
   do i = 1,NLONS
      do j = 1, NLATS
         LS%TAIR(i,j,1,NTIME) = t_ref2m((j-1)*NLONS+i)
      enddo
   enddo

   print*, '2m air temperature is',LS%TAIR(1,1,1,NTIME),LS%TAIR(NLONS,NLATS,1,NTIME)
   
!--------------------------------------------------------------------

   print*, 'CMEM input read Successfully close!'
   
   ! --------------------------------------------
   ! Declaring global variables for CMEM which
   ! are needed for storing final variables as
   ! NetCDF ouput:
   NLONS_SM = NLONS
   NLATS_SM = NLATS
   NLVLS_SM = 1_JPIM
   NTIMES_SM = NTIME
   IF( ALLOCATED (xlons)) DEALLOCATE (xlons)
   ALLOCATE(xlons(NLONS))
   xlons = LS%lons
   IF( ALLOCATED (xlats)) DEALLOCATE (xlats)
   ALLOCATE(xlats(NLATS))
   xlats = LS%lats
   IF( ALLOCATED (xlvls)) DEALLOCATE (xlvls)
   ALLOCATE(xlvls(nlev))
   xlvls = LS%levsoi
   IF( ALLOCATED (xtimes)) DEALLOCATE (xtimes)
   ALLOCATE(xtimes(NTIME))
   xtimes = LS%time
   
   deallocate(idxlev)
   return
!#endif
  
end subroutine read_CLM_pdaf
! ===================== END OF READING CLM FILE SUBROUTINE =========

end module rdclm_wrcmempdaf

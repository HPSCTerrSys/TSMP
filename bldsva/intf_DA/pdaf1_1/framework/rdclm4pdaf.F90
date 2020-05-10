module rdclm4pdaf 

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
  
  subroutine read_CLM_pdaf(LS,SAT)

!#if model == tag_model_clm
   USE PARKIND1,      ONLY : JPIM
   USE clmtype      , only : clm3,nameg,namec,namep
   USE clm_varpar   , only : nlevsoi,nlevsno!,nlevgrnd
   use shr_kind_mod,  only: r8 => shr_kind_r8
   use clm4cmem,      only: SATELLITE, local_incidence, gradientm
   use domainMod,     only: ldomain
   use decompMod,     only: get_proc_global,get_proc_bounds,adecomp
   ! use surfrdMod, only: pctpft
   use mod_parallel_pdaf, only: COMM_filter
   ! use mod_parallel_model, only: MPI_COMM_WORLD
   use mod_tsmp,      only: idx_map_subvec2state_fortran,pf_statevecsize
   ! use clm_varsur,    only: longxy, latixy
   ! use mod_assimilation,  ONLY:longxy, latixy! obs_index_p
   USE YOMCMEMPAR,    ONLY: LGPRINT
   USE YOMLUN,        ONLY: NULOUT
   !USE oas_clm_vardef
   USE spmdMod      , only : masterproc,iam,mpicom,npes
   USE spmdGathScatMod , only : gather_data_to_master
   implicit none
   type(SATELLITE), intent(inout) :: SAT
   type(CLM_DATA), intent(out) :: LS
   integer(kind=JPIM) :: Nvars(3)
   integer, dimension(:), allocatable :: idxlev
   
   INTEGER :: numg           ! total number of gridcells across all processors
   INTEGER :: numl           ! total number of landunits across all processors
   INTEGER :: numc           ! total number of columns across all processors
   INTEGER :: nump           ! total number of pfts across all processors
   integer :: begp, endp   ! per-proc beginning and ending pft indices
   integer :: begc, endc   ! per-proc beginning and ending column indices
   integer :: begl, endl   ! per-proc beginning and ending landunit indices
   integer :: begg, endg   ! per-proc gridcell ending gridcell indices
   integer :: clm_begg,clm_endg 
   integer :: nlev, idxtime=1
   
   integer :: NINC, NLATS, NLONS, NTIME, NLEVPFT, &
             ier,gi,i,j,k_lev,cc=1,offset=0,nlevgrnd

   real(r8), pointer :: latdeg(:)       !latitude (degrees)
   real(r8), pointer :: londeg(:)       !longitude (degrees)
   !real(r8), pointer :: dz(:,:)         !layer thickness (m)  (-nlevsno+1:nlevgrnd) 
   real(r8), pointer :: z(:,:)          !layer depth (m) (-nlevsno+1:nlevgrnd) 
   real(r8), pointer :: swc(:,:)
   real(r8), pointer :: psand(:,:)
   real(r8), pointer :: pclay(:,:)
   real(r8), pointer :: snowdp(:)
   real(r8), pointer :: snowliq(:)
   real(r8), pointer :: frac_sno(:)
   real(r8), pointer :: t_soisno(:,:)
   real(r8), pointer :: t_grnd(:)
   real(r8), pointer :: tlai(:)        
   real(r8), pointer :: t_ref2m(:)
   integer , pointer :: itype(:)
   real(r8), pointer :: topo(:)
   real(r8), pointer :: time(:)
   real(r8), pointer :: fcov(:)

   real(r8), pointer :: latdeg_PE(:)       !latitude (degrees)
   real(r8), pointer :: londeg_PE(:)       !longitude (degrees)
   real(r8), pointer :: z_PE(:,:)          !layer depth (m) (-nlevsno+1:nlevgrnd) 
   real(r8), pointer :: swc_PE(:,:)
   real(r8), pointer :: psand_PE(:,:)
   real(r8), pointer :: pclay_PE(:,:)
   real(r8), pointer :: snowdp_PE(:)
   real(r8), pointer :: snowliq_PE(:)
   real(r8), pointer :: frac_sno_PE(:)
   real(r8), pointer :: t_soisno_PE(:,:)
   real(r8), pointer :: t_grnd_PE(:)
   real(r8), pointer :: tlai_PE(:)
   real(r8), pointer :: t_ref2m_PE(:)
   integer , pointer :: itype_PE(:)
   real(r8), pointer :: topo_PE(:)
   !real(r8), pointer :: time(:)
   real(r8), pointer :: fcov_PE(:) 
   real(r8), pointer :: latdeg_tmp(:)       !latitude (degrees)
   real(r8), pointer :: londeg_tmp(:)       !longitude (degrees)
   real(r8), pointer :: z_tmp(:,:)          !layer depth (m) (-nlevsno+1:nlevgrnd) 
   real(r8), pointer :: swc_tmp(:,:)
   real(r8), pointer :: psand_tmp(:,:)
   real(r8), pointer :: pclay_tmp(:,:)
   real(r8), pointer :: snowdp_tmp(:)
   real(r8), pointer :: snowliq_tmp(:)
   real(r8), pointer :: frac_sno_tmp(:)
   real(r8), pointer :: t_soisno_tmp(:,:)
   real(r8), pointer :: t_grnd_tmp(:)
   real(r8), pointer :: tlai_tmp(:)
   real(r8), pointer :: t_ref2m_tmp(:)
   integer , pointer :: itype_tmp(:)
   real(r8), pointer :: topo_tmp(:)
   ! real(r8), pointer :: time(:)
   real(r8), pointer :: fcov_tmp(:)
   
   latdeg_PE   => clm3%g%latdeg 
   londeg_PE   => clm3%g%londeg 
   z_PE        => clm3%g%l%c%cps%z
   swc_PE      => clm3%g%l%c%cws%h2osoi_vol
   psand_PE    => clm3%g%l%c%cps%psand
   pclay_PE    => clm3%g%l%c%cps%pclay
   snowdp_PE   => clm3%g%l%c%cps%snowdp
   snowliq_PE  => clm3%g%l%c%cws%snowliq
   frac_sno_PE => clm3%g%l%c%cps%frac_sno
   t_soisno_PE => clm3%g%l%c%ces%t_soisno
   t_grnd_PE   => clm3%g%l%c%ces%t_grnd
   tlai_PE     => clm3%g%l%c%p%pps%tlai
   t_ref2m_PE  => clm3%g%l%c%p%pes%t_ref2m
   itype_PE    => clm3%g%l%c%p%itype 
   time        => clm3%g%l%c%p%pepv%days_active
   fcov_PE     => clm3%g%l%c%cws%fcov         
   topo_PE     => ldomain%topo
   
   CALL get_proc_global(numg,numl,numc,nump)
   call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

   IF (LGPRINT) WRITE(NULOUT,*)  "begg/endg is ",begg,endg
   clm_begg     = 1
   clm_endg     = adecomp%gdc2i(numg)*adecomp%gdc2j(numg)
   NLATS        = adecomp%gdc2j(numg)
   NLONS        = adecomp%gdc2i(numg)
   NINC         = 1 
   NTIME        = 1
   NLEVPFT      = 1 
   nlevgrnd     = nlevsoi+nlevsno
   IF (LGPRINT) WRITE(NULOUT,*)  "clm_begg/clm_endg is",clm_begg,clm_endg,masterproc
  
   IF (masterproc) THEN
   !   IF (.NOT. associated(swc)) THEN
         ALLOCATE(latdeg(numg), stat=ier)
         ALLOCATE(londeg(numg), stat=ier)
         ALLOCATE(z(nlevgrnd,numg), stat=ier)
         ALLOCATE(swc(nlevsoi,numg), stat=ier)
         ALLOCATE(psand(nlevsoi,numg), stat=ier)
         ALLOCATE(pclay(nlevsoi,numg), stat=ier)
         ALLOCATE(snowdp(numg), stat=ier)
         ALLOCATE(snowliq(numg), stat=ier)
         ALLOCATE(frac_sno(numg), stat=ier)
         ALLOCATE(t_soisno(nlevgrnd,numg), stat=ier)
         ALLOCATE(t_grnd(numg), stat=ier)
         ALLOCATE(tlai(numg), stat=ier)
         ALLOCATE(t_ref2m(numg), stat=ier)
         ALLOCATE(itype(numg), stat=ier)
         ALLOCATE(topo(numg), stat=ier)
         ALLOCATE(fcov(numg), stat=ier)   
   !   END IF
   END IF
 
   ALLOCATE(latdeg_tmp(begg:endg), stat=ier)
   ALLOCATE(londeg_tmp(begg:endg), stat=ier)
   ALLOCATE(z_tmp(nlevgrnd,begg:endg), stat=ier)
   ALLOCATE(swc_tmp(nlevsoi,begg:endg), stat=ier)
   ALLOCATE(psand_tmp(nlevsoi,begg:endg), stat=ier)
   ALLOCATE(pclay_tmp(nlevsoi,begg:endg), stat=ier)
   ALLOCATE(snowdp_tmp(begg:endg), stat=ier)
   ALLOCATE(snowliq_tmp(begg:endg), stat=ier)
   ALLOCATE(frac_sno_tmp(begg:endg), stat=ier)
   ALLOCATE(t_soisno_tmp(nlevgrnd,begg:endg), stat=ier)
   ALLOCATE(t_grnd_tmp(begg:endg), stat=ier)
   ALLOCATE(tlai_tmp(begg:endg), stat=ier)
   ALLOCATE(t_ref2m_tmp(begg:endg), stat=ier)
   ALLOCATE(itype_tmp(begg:endg), stat=ier)
   ALLOCATE(topo_tmp(begg:endg), stat=ier)
   ALLOCATE(fcov_tmp(begg:endg), stat=ier)
 
   !latdeg_tmp = TRANSPOSE(latdeg_PE)
   latdeg_tmp   = latdeg_PE
   CALL gather_data_to_master(latdeg_tmp, latdeg, clmlevel=nameg)
   londeg_tmp   = londeg_PE 
   CALL gather_data_to_master(londeg_tmp, londeg, clmlevel=nameg)
   !dz_tmp       = TRANSPOSE(dz_PE)
   !CALL gather_data_to_master(dz_tmp, dz, clmlevel=namec)
   z_tmp        = TRANSPOSE(z_PE)
   CALL gather_data_to_master(z_tmp, z, clmlevel=namec)
   swc_tmp      = TRANSPOSE(swc_PE)
   CALL gather_data_to_master(swc_tmp, swc, clmlevel=namec)
   psand_tmp    = TRANSPOSE(psand_PE)
   CALL gather_data_to_master(psand_tmp, psand, clmlevel=namec)
   pclay_tmp    = TRANSPOSE(pclay_PE)
   CALL gather_data_to_master(pclay_tmp, pclay, clmlevel=namec)
   snowdp_tmp   = snowdp_PE
   CALL gather_data_to_master(snowdp_tmp,snowdp, clmlevel=namec)
   snowliq_tmp  = snowliq_PE
   CALL gather_data_to_master(snowliq_tmp, snowliq, clmlevel=namec)
   frac_sno_tmp   = frac_sno_PE
   CALL gather_data_to_master(frac_sno_tmp,frac_sno, clmlevel=namec)
   t_soisno_tmp = TRANSPOSE(t_soisno_PE)
   CALL gather_data_to_master(t_soisno_tmp, t_soisno, clmlevel=namec)
   t_grnd_tmp   = t_grnd_PE
   CALL gather_data_to_master(t_grnd_tmp, t_grnd, clmlevel=namec)
   tlai_tmp     = tlai_PE
   CALL gather_data_to_master(tlai_tmp, tlai, clmlevel=namep)
   t_ref2m_tmp  = t_ref2m_PE
   CALL gather_data_to_master(t_ref2m_tmp, t_ref2m, clmlevel=namep)
   itype_tmp    = itype_PE
   CALL gather_data_to_master(itype_tmp, itype, clmlevel=namec)
   fcov_tmp     = fcov_PE
   CALL gather_data_to_master(fcov_tmp, fcov, clmlevel=namec)
   topo_tmp     = topo_PE
   CALL gather_data_to_master(topo_tmp, topo, clmlevel=nameg)
   
   !write(*,*) 'whether it is the masterproc', masterproc
   !write(*,*) 'kl_comm are',kl_comm,iam
   !CALL MPI_Barrier(MPI_COMM_WORLD, nerror)
   IF (LGPRINT) WRITE(NULOUT,*) 'mpicom/iam/npes are',mpicom,iam,npes
   CALL MPI_Barrier(mpicom, ier)   
 
   IF (masterproc) THEN
   ! Get the varids of the latitude and longitude coordinate variables.
   ! Read the latitude and longitude data.
   allocate(LS%lats(NLATS))
   allocate(LS%lons(NLONS))
   do gi = 1,numg
      i = adecomp%gdc2i(gi)
      j = adecomp%gdc2j(gi)
      if(i==1)   LS%lats(j) = latdeg(gi)
      if(j==1)   LS%lons(i) = londeg(gi)
   end do 
   IF (LGPRINT) WRITE(NULOUT,*) 'input LS LAT/LON begins at ',LS%lats(1),LS%lons(1)
   IF (LGPRINT) WRITE(NULOUT,*) 'input LS LAT/LON ends at ',LS%lats(NLATS),LS%lons(NLONS)
   !WRITE(NULOUT,*) LS%lats,LS%lons   

   ! Setting default values for input variables in case not passed in call  
   nlev         = nlevsoi
   IF (LGPRINT) WRITE(NULOUT,*) 'input LS NLEV is ',nlev
   allocate(idxlev(nlev))
   idxlev = (/(I,I=1,nlev)/)
   if(.not.allocated(LS%levsoi)) allocate(LS%levsoi(nlev))
   DO k_lev = 1,nlev
         LS%levsoi(k_lev) = z(k_lev+nlevgrnd-nlevsoi,1)
   END DO
   IF (LGPRINT) WRITE(NULOUT,*) 'LSM layers depth are '
   IF (LGPRINT) WRITE(NULOUT,'(5f9.3)') LS%levsoi

   LS%DATESTRING = '18/04/20 10:28:22'
   if(.not.allocated(LS%time)) allocate(LS%time(NTIME))
   LS%time = time(idxtime)
   NINC = size(SAT%theta)
   if(allocated(SAT%time)) deallocate(SAT%time)
   allocate(SAT%time(NTIME))
   SAT%time = time(idxtime)
   IF (LGPRINT) WRITE(NULOUT,*) 'SAT time is ',SAT%time 
   
   ! 0 -> Reading TOPO:
   allocate(LS%Z(NLONS,NLATS))
   ! geopotential height from m(CLM) to km(CMEM)
   do gi = 1,numg
      i = adecomp%gdc2i(gi)
      j = adecomp%gdc2j(gi)
      LS%Z(i,j) = topo(gi)/1000
   enddo
   IF (LGPRINT) WRITE(NULOUT,*) 'geopotential height(km) is ',LS%Z(1,1),LS%Z(NLONS,NLATS)

   ! 0.1 -> reading the longitude and latitude grid-boxes:
   allocate(LS%latixy(NLONS,NLATS), LS%longxy(NLONS,NLATS))
   !LS%latixy = reshape(latixy, (/NLONS,NLATS/))
   !LS%longxy = reshape(longxy, (/NLONS,NLATS/))
   do gi = 1,numg
      i = adecomp%gdc2i(gi)
      j = adecomp%gdc2j(gi)
      LS%latixy(i,j) = latdeg(gi)
      LS%longxy(i,j) = londeg(gi)
   enddo

   IF (LGPRINT) WRITE(NULOUT,*) 'reading the longitude and latitude is done'
   IF (LGPRINT) WRITE(NULOUT,*) LS%latixy(1,1),LS%latixy(NLONS,NLATS)   
   IF (LGPRINT) WRITE(NULOUT,*) LS%longxy(1,1),LS%longxy(NLONS,NLATS)
   ! 0.2 -> calculating slope and aspect from TOPO Z:
   allocate(LS%slope(NLONS,NLATS), LS%aspect(NLONS,NLATS))
   allocate(LS%theta_inc(NLONS,NLATS,NINC))
   call gradientm(LS%Z, LS%longxy, LS%latixy, LS%slope,&
        LS%aspect, LS%resol_km)

   ! 0.3 -> calculating pixel-wise incidence angle theta_inc:
   LS%theta = SAT%theta
   LS%theta_inc = local_incidence(LS%slope, LS%aspect, SAT)
   IF (LGPRINT) WRITE(NULOUT,*) 'theta is ',LS%theta
   
   ! 1. Making up soil moisture [%/%]        
   allocate(LS%SWVL(NLONS,NLATS,nlev,NTIME))
   !LS%SWVL = reshape(swc, (/NLONS,NLATS,NTIME,nlev/))
   do gi = 1,numg
      i = adecomp%gdc2i(gi)
      j = adecomp%gdc2j(gi)
      do k_lev = 1,nlev
         if (swc(k_lev,gi).EQ.1E+36) then
            LS%SWVL(i,j,k_lev,NTIME) = -1E+34
         else
            LS%SWVL(i,j,k_lev,NTIME) = swc(k_lev,gi)
         end if 
      end do
   end do
   IF (LGPRINT) WRITE(NULOUT,*) 'soil moisture is ',LS%SWVL(1,1,1,1),LS%SWVL(NLONS,NLATS,nlev,NTIME)
   
   ! 2. Making up soil temperature [K]
   allocate(LS%STL(NLONS,NLATS,nlev,NTIME))  
   !LS%STL = reshape(t_soisno,(/NLONS,NLATS,NTIME,nlev/))
   do gi = 1,numg
      i = adecomp%gdc2i(gi)
      j = adecomp%gdc2j(gi)
      do k_lev = 1,nlev
         if (t_soisno(k_lev+nlevgrnd-nlevsoi,gi).EQ.1E+36) then
            LS%STL(i,j,k_lev,NTIME) = -1E+34
         else
            LS%STL(i,j,k_lev,NTIME) = t_soisno(k_lev+nlevgrnd-nlevsoi,gi)
         end if
      end do
   end do
   IF (LGPRINT) WRITE(NULOUT,*) 'soil temperature is',LS%STL(1,1,1,1),LS%STL(NLONS,NLATS,nlev,NTIME)
   
   ! 3. Making up snow depth
   allocate(LS%SD(NLONS,NLATS,1,NTIME))
   !LS%SD = reshape(snowdp, (/NLONS,NLATS,1,NTIME/))
   do gi = 1,numg
      i = adecomp%gdc2i(gi)
      j = adecomp%gdc2j(gi)
      if (snowdp(gi).EQ.1E+36) then
         LS%SD(i,j,1,NTIME) = -1E+34
      else
         LS%SD(i,j,1,NTIME) = snowdp(gi)
      end if
   end do 
   IF (LGPRINT) WRITE(NULOUT,*) 'snow depth is ',LS%SD(1,1,1,NTIME),LS%SD(NLONS,NLATS,1,NTIME)
   
   ! 4. Making up snow density [kg/m^3]
   !  water equivalent of snow/snow depth [m]
   allocate(LS%RSN(NLONS,NLATS,1,NTIME))
   do gi = 1,numg
      i = adecomp%gdc2i(gi)
      j = adecomp%gdc2j(gi)
      if (snowdp(gi).EQ.1E+36.OR.snowdp(gi).LE.0.0) then
         LS%RSN(i,j,1,NTIME) = -1E+34
      elseif (snowliq(gi).LE.0.0.OR.snowliq(gi).EQ.1E+36) then
         LS%RSN(i,j,1,NTIME) = -1E+34
      else
         LS%RSN(i,j,1,NTIME) = snowliq(gi)/snowdp(gi)
      end if
      if (LS%RSN(i,j,1,NTIME).LE.0.0) then
         LS%RSN(i,j,1,NTIME) = 0.0
      end if
   end do
   IF (LGPRINT) WRITE(NULOUT,*) 'snow RSN is ',LS%RSN(1,1,1,NTIME),LS%RSN(NLONS,NLATS,1,NTIME)
   
   ! 5. Making up skin temperature [K]
   allocate(LS%TSKIN(NLONS,NLATS,1,NTIME))
   ! LS%TSKIN = reshape(t_grnd, (/NLONS,NLATS,1,NTIME/))
   ! IF (LGPRINT) WRITE(NULOUT,*) t_grnd
   do gi = 1,numg
      i = adecomp%gdc2i(gi)
      j = adecomp%gdc2j(gi)
      if (t_grnd(gi).EQ.1E+36) then
         LS%TSKIN(i,j,1,NTIME) = -1E+34
      else
         LS%TSKIN(i,j,1,NTIME) = t_grnd(gi)
      end if
   end do 
   IF (LGPRINT) WRITE(NULOUT,*) 'skin temperature is',LS%TSKIN(1,1,1,NTIME),LS%TSKIN(NLONS,NLATS,1,NTIME) 
   
   ! 6 texture
   allocate(LS%CLAY(NLONS,NLATS))
   allocate(LS%SAND(NLONS,NLATS))
   do gi = 1,numg
      i = adecomp%gdc2i(gi)
      j = adecomp%gdc2j(gi)
      LS%SAND(i,j) = psand(1,gi)
      LS%CLAY(i,j) = pclay(1,gi)
   enddo
   IF (LGPRINT) WRITE(NULOUT,*) 'sand fraction is ',LS%SAND(1,1),LS%SAND(NLONS,NLATS)
   IF (LGPRINT) WRITE(NULOUT,*) 'clay fraction is ',LS%CLAY(1,1),LS%CLAY(NLONS,NLATS)
   ! 7. vegetation
 
   ! Vegetation
   ! Vegetation types    | CLM VR0  | ECOCLIMAP | TESSEL :
   ! ---------------------------------------------------------------
   ! bare_soil           | 0        |           |
   ! needle_leaf_forest  | 1        |  2        | 3      (High Veg)
   ! broad_leaf_forest   | 7        |  1        | 4      (High Veg)
   ! grassland           | 14       |  5        | 2      (Low Veg)
   ! crops/agricul_land  | 15       |  6        | 1      (Low Veg)
   ! ---------------------------------------------------------------
   ! -------------------LSN new PFT table---------------------------
   ! PFT types                           | CLM3.5      | TESSEL
   ! ---------------------------------------------------------------
   ! not vegetated                       | 0           | 
   ! needleleaf evergreen temperate tree | 1           | 3   (High Veg)
   ! needleleaf evergreen boreal tree    | 2           | 3   (High Veg)
   ! needleleaf deciduous boreal tree    | 3           | 4   (High Veg)
   ! broadleaf evergreen tropical tree   | 4           | 6   (High Veg)
   ! broadleaf evergreen temperate tree  | 5           | 6   (High Veg)
   ! broadleaf deciduous tropical tree   | 6           | 5   (High Veg)
   ! broadleaf deciduous temperate tree  | 7           | 5   (High Veg)
   ! broadleaf deciduous boreal tree     | 8           | 5   (High Veg)
   ! broadleaf evergreen shrub           | 9           | 16  (Low Veg)
   ! broadleaf deciduous temperate shrub | 10          | 17  (Low Veg)
   ! broadleaf deciduous boreal shrub    | 11          | 17  (Low Veg)
   ! c3 arctic grass                     | 12          | 2   (Low Veg)
   ! c3 non-arctic grass                 | 13          | 2   (Low Veg)
   ! c4 grass                            | 14          | 7   (Low Veg)
   ! corn                                | 15          | 1   (Low Veg)
   ! wheat                               | 16          | 1   (Low Veg)
   ! ---------------------------------------------------------------
     
   ! 7.1 ! 7.1 -> high/low vegetation fraction. TerrsysMP is high spatial
   ! resolution, so one pft for each cell. Set fraction to 1.  
  
   allocate(LS%CVL(NLONS,NLATS))   ! TESSEL fraction of low veg
   allocate(LS%CVH(NLONS,NLATS))   ! TESSEL fraction of high veg
    
   do gi = 1,numg
      i = adecomp%gdc2i(gi)
      j = adecomp%gdc2j(gi)
      ! High veg
      if(itype(gi).ge.1 .and. itype(gi).le.8) then
      !      LS%CVH(i,j) = (1.0-fcov(gi))*frac_sno(gi)
            LS%CVH(i,j) = 1.0
            LS%CVL(i,j) = 0
      ! Low veg
      elseif(itype(gi).ge.9) then
            LS%CVH(i,j) = 0
            LS%CVL(i,j) = 1.0
      !      LS%CVL(i,j) = (1.0-fcov(gi))*frac_sno(gi)
      ! Bare soil
      else
            LS%CVH(i,j) = 0
            LS%CVL(i,j) = 0
      endif
   end do
    
   ! 7.2 high/low vegetation type 
   allocate(LS%TVL(NLONS,NLATS))   ! type of low veg
   allocate(LS%TVH(NLONS,NLATS))   ! type of high veg 
   do gi = 1,numg
      i = adecomp%gdc2i(gi)
      j = adecomp%gdc2j(gi)
      ! High veg
      if(itype(gi).eq.1.or.itype(gi).eq.2) then
         LS%TVH(i,j)   = 3
         LS%TVL(i,j)   = 0
      elseif(itype(gi).eq.3) then
         LS%TVH(i,j)   = 4
         LS%TVL(i,j)   = 0
      elseif(itype(gi).eq.4.or.itype(gi).eq.5) then
         LS%TVH(i,j)   = 6
         LS%TVL(i,j)   = 0
      elseif(itype(gi).eq.6.or.itype(gi).eq.7.or.itype(gi).eq.8) then
         LS%TVH(i,j)   = 5
         LS%TVL(i,j)   = 0
      ! Low veg
      elseif(itype(gi).eq.9) then
         LS%TVH(i,j)   = 0
         LS%TVL(i,j)   = 16
      elseif(itype(gi).eq.10.or.itype(gi).eq.11) then
         LS%TVH(i,j)   = 0
         LS%TVL(i,j)   = 17
      elseif(itype(gi).eq.12.or.itype(gi).eq.13) then
         LS%TVH(i,j)   = 0
         LS%TVL(i,j)   = 2
      elseif(itype(gi).eq.14) then
         LS%TVH(i,j)   = 0
         LS%TVL(i,j)   = 7
      elseif(itype(gi).eq.15.or.itype(gi).eq.16) then
         LS%TVH(i,j)   = 0
         LS%TVL(i,j)   = 1
      ! Bare soil
      else
         LS%TVH(i,j)   = 0
         LS%TVL(i,j)   = 0
      endif
   end do

   ! 7.3 water fraction [-] / land franction [-]
   if(allocated(LS%WATER)) deallocate(LS%WATER)
   allocate(LS%WATER(NLONS,NLATS))
   if(allocated(LS%LSM)) deallocate(LS%LSM)
   allocate(LS%LSM(NLONS,NLATS)) 
   do gi = 1,numg
      i = adecomp%gdc2i(gi)
      j = adecomp%gdc2j(gi)
      !if(fcov(gi).gt.0) then
         LS%WATER(i,j) = fcov(gi)
         LS%LSM(i,j)   = 1.0-fcov(gi)
      !else
      !   LS%WATER(i,j) = 0
      !   LS%LSM(i,j)   = 1
      !endif
   end do

   ! 7.4 LAI of low vegetation
   if(allocated(LS%LAIL)) deallocate(LS%LAIL)
   allocate(LS%LAIL(NLONS,NLATS,1,NTIME))
   do gi = 1,numg
      i = adecomp%gdc2i(gi)
      j = adecomp%gdc2j(gi)
      LS%LAIL(i,j,1,NTIME) = tlai(gi)
   enddo
   IF (LGPRINT) WRITE(NULOUT,*) 'LAI of low vegeation is ',LS%LAIL(1,1,1,NTIME),LS%LAIL(NLONS,NLATS,1,NTIME)

   ! 8. 2m Air Temperature [K]
   if (allocated(LS%TAIR)) deallocate(LS%TAIR)
   allocate(LS%TAIR(NLONS,NLATS,1,NTIME))
   do gi = 1,numg
      i = adecomp%gdc2i(gi)
      j = adecomp%gdc2j(gi)
      if (t_ref2m(gi).EQ.1E+36) then
         LS%TAIR(i,j,1,NTIME) = -1E+34
      else
         LS%TAIR(i,j,1,NTIME) = t_ref2m(gi)
      end if
   enddo

   IF (LGPRINT) WRITE(NULOUT,*) '2m air temperature is',LS%TAIR(1,1,1,NTIME),LS%TAIR(NLONS,NLATS,1,NTIME)
   
   IF (LGPRINT) WRITE(NULOUT,*) 'CMEM input read Successfully close!'
   
   ! Declaring global variables for CMEM which are needed for storing final variables as
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
   DEALLOCATE(latdeg,londeg,z,swc,psand,pclay,snowdp,snowliq,t_soisno,t_grnd,tlai,t_ref2m,itype,topo,fcov)

END IF ! if masterproc
  
DEALLOCATE(latdeg_tmp,londeg_tmp,z_tmp,swc_tmp,psand_tmp,pclay_tmp,snowdp_tmp,snowliq_tmp)
DEALLOCATE(t_soisno_tmp,t_grnd_tmp,tlai_tmp,t_ref2m_tmp,itype_tmp,topo_tmp,fcov_tmp)

end subroutine read_CLM_pdaf
! ===================== END OF READING CLM FILE SUBROUTINE =========

end module rdclm4pdaf

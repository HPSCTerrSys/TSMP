SUBROUTINE send_fld_2pfl

!---------------------------------------------------------------------------------
! Description:
! This subroutine sends the fluxes from CLM3.5 to ParFlow
!---------------------------------------------------------------------------------
! Current Code Owner: TR32, Z4: Mauro Sulis
!    phone: 02287360704
!    email: msulis@uni-bonn.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        2012/01/30 Mauro Sulis
!   Modfied and Implemented in CLM3.5, Initial release
! 1.1        2013/06/30 Mauro Sulis
! Bug fix in dz for sending fluxes with multiple threads
! @VERSION@    @DATE@     <Your name>
!  <Modification comments>         
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
! ============================================================================

! Declarations:
!
! Modules used:

 
USE shr_kind_mod , only : r8 => shr_kind_r8
USE clm_atmlnd   , only : clm_l2a, atm_l2a, clm_mapl2a
USE clmtype      , only : clm3, nameg
USE subgridAveMod, only : p2g, c2g
USE domainMod    , only : latlon_type
USE clm_varpar   , only : nlevsoi
USE decompMod    , only : get_proc_global, get_proc_bounds, adecomp
USE spmdGathScatMod , only : gather_data_to_master
USE spmdMod      , only : masterproc
USE clm_time_manager        , only :   &
                                        get_nstep,        &! return timestep number
                                        dtime, nelapse              ! timestep in second
USE oas_clm_vardef
USE netcdf

! ============================================================================

IMPLICIT NONE

! ============================================================================
   
INTEGER :: numg           ! total number of gridcells across all processors
INTEGER :: numl           ! total number of landunits across all processors
INTEGER :: numc           ! total number of columns across all processors
INTEGER :: nump           ! total number of pfts across all processors
INTEGER :: begg,endg      ! local beg/end gridcells gdc
INTEGER :: begl,endl      ! local beg/end landunits
INTEGER :: begc,endc      ! local beg/end columns
INTEGER :: begp,endp      ! local beg/end pfts

INTEGER ::   isec, info, jn, jj, ji, g1, jx, i    ! temporary integer

REAL(r8), POINTER :: dz(:,:)           !layer depth (m)
REAL(r8), POINTER :: dz_tmp(:,:)       !layer depth (m)
REAL(r8), POINTER :: dz_snd(:,:)      !layer depth (m)
   
REAL(r8), ALLOCATABLE :: snd_field(:,:)
REAL(r8), POINTER :: snd_tmp(:,:)

REAL(r8), ALLOCATABLE ::   fsnd(:,:,:)            ! temporary arrays

REAL(r8), POINTER :: pfl_flx_total_col(:,:)
REAL(r8), POINTER :: pfl_flx_total_gcell(:,:)
REAL(r8), POINTER :: pfl_flx_total_gcell_tmp(:,:)
!
INTEGER ,DIMENSION(4) :: dimids
INTEGER ,DIMENSION(1) :: il_var_id
INTEGER :: il_file_id, ncvarid(1), status

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine send_fld2pfl 
!------------------------------------------------------------------------------

! Set pointer to dz
 dz            => clm3%g%l%c%cps%dz

! Set pointers into derived type (column level)
   
 pfl_flx_total_col    => clm3%g%l%c%cwf%pfl_flx_total

! Set pointers into derived type (gridcell level)

pfl_flx_total_gcell   => clm3%g%gwf%pfl_flx_total_gcell


! Get total global number of grid cells, landunits, columns and pfts

CALL get_proc_global(numg,numl,numc,nump)
CALL get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

! allocate global clm buffers
IF (masterproc) THEN
   IF (.NOT. associated(snd_tmp)) THEN
      ALLOCATE(snd_tmp(nlevsoi,numg), stat=nerror)
      IF (nerror /= 0) THEN
         CALL prism_abort_proto( ncomp_id, 'send_fld_2pf', 'Failure in allocating snd_tmp' )
         RETURN
      END IF
      ALLOCATE(snd_field(nlevsoi,numg), stat=nerror)
      IF (nerror /= 0) then
      CALL prism_abort_proto( ncomp_id, 'send_fld_2pf', 'Failure in allocating snd_field' )
      RETURN
      END IF
      ALLOCATE(dz_snd(nlevsoi,numg), stat=nerror)
      IF (nerror /= 0) then
      CALL prism_abort_proto( ncomp_id, 'send_fld_2pf', 'Failure in allocating dz_snd' )
      RETURN
      END IF
   END IF
END IF

ALLOCATE(pfl_flx_total_gcell_tmp(nlevsoi,begg:endg), stat=nerror) !CMS: We make the allocation at local level
IF (nerror /= 0) THEN
    CALL prism_abort_proto( ncomp_id, 'send_fld_2pf', 'Failure in allocating snd_tmp' )
    RETURN
END IF
   
ALLOCATE(dz_tmp(nlevsoi,begg:endg), stat=nerror)
IF (nerror /= 0) then
CALL prism_abort_proto( ncomp_id, 'send_fld_2pf', 'Failure in allocating dz_tmp' )
RETURN
END IF
   
CALL c2g(begc, endc, begl, endl, begg, endg, nlevsoi, pfl_flx_total_col, pfl_flx_total_gcell, &
         c2l_scale_type= 'unity', l2g_scale_type='unity')
  
   
pfl_flx_total_gcell_tmp = TRANSPOSE(pfl_flx_total_gcell) !CMS: We send the transpose of the array in order to match
                                                         !the 2d dimensions in subroutine gather_2darray_real

dz_tmp = TRANSPOSE(dz(:,1:nlevsoi))                      !CMS: dz (:,1:nlevsoi) as it considers snow levels

CALL gather_data_to_master(pfl_flx_total_gcell_tmp, snd_tmp, clmlevel=nameg)

CALL gather_data_to_master(dz_tmp, dz_snd, clmlevel=nameg)
   
IF (masterproc) snd_field=snd_tmp

IF (masterproc) THEN

   ji = adecomp%gdc2i(numg)
   jj = adecomp%gdc2j(numg)
      
ALLOCATE ( fsnd(ndlon,ndlat,nlevsoi), stat=nerror)
   IF (nerror /= 0) THEN
      CALL prism_abort_proto( ncomp_id, 'send_fld_2pf', 'Failure in allocating fsnd' )
      RETURN
   END IF

! zero on unmasked points
   fsnd = 0._r8

   DO jn = 1, vsnd
   DO g1 = 1, numg
        ji = adecomp%gdc2i(g1)
        jj = adecomp%gdc2j(g1)
        fsnd(ji,jj,jn) = snd_field(jn,g1)*3.6_r8/dz_snd(jn,g1)!CMS: Note that we got the transpose from gather_data_to master
                                                              !     We assume timing for pfl is hours and seconds for clm, 
                                                              !     length is m for pfl and mm for clm
   END DO
   END DO

   isec = dtime * ( get_nstep() -1 )
   
!------------------------------------------------------------------------------
   IF ( IOASISDEBUGLVL == 1 ) THEN
   
   IF (isec == 0) THEN
    status =  nf90_create("debugsnd_clm_pfl.nc", NF90_CLOBBER, il_file_id)
    status =  nf90_def_dim(il_file_id, "x", ndlon, dimids(1))
    status =  nf90_def_dim(il_file_id, "y", ndlat, dimids(2))
    status =  nf90_def_dim(il_file_id, "z", nlevsoi, dimids(3))
    status =  nf90_def_dim(il_file_id, "t", nelapse, dimids(4))
    status =  nf90_def_var(il_file_id, "CLMFLX", NF90_DOUBLE, dimids, ncvarid(1))
    status =  nf90_enddef(il_file_id)
    status =  nf90_close(il_file_id)
   END IF
    
    status = nf90_open("debugsnd_clm_pfl.nc", NF90_WRITE, il_file_id)
    status = nf90_inq_varid(il_file_id, "CLMFLX" , ncvarid(1))
    status = nf90_put_var( il_file_id, ncvarid(1), fsnd(:,:,:),  &
                          start = (/ 1, 1, 1, get_nstep()/),          &
                          count = (/ ndlon, ndlat, nlevsoi, 1 /) )
    status = nf90_close(il_file_id)
   ENDIF           
!------------------------------------------------------------------------------

   DO jn = 1, vsnd

   IF( ssnd(100+jn)%laction )  THEN !CMS the first 100th snd fields are allocated for COSMO-CLM coupling
    

   CALL oas_clm_snd( 100+jn, isec, fsnd(:,:,jn), info ) !CMS the first 100th snd fields are allocated for COSMO-CLM coupling

    
   END IF !if ssnd is TRUE

   END DO ! jn 
   
   isec = dtime * get_nstep()

   DEALLOCATE(fsnd)
   DEALLOCATE(snd_field)
   DEALLOCATE(snd_tmp)
   DEALLOCATE(dz_snd)

   END IF ! if masterproc

DEALLOCATE(pfl_flx_total_gcell_tmp)
DEALLOCATE(dz_tmp)

CALL MPI_Barrier(kl_comm, nerror)

!CPS WRITE(6,*) 'oas_snd_pfl complete ...'

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE send_fld_2pfl

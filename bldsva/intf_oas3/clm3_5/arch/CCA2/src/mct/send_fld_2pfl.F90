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
USE spmdMod      , only : masterproc, iam
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

REAL(r8), ALLOCATABLE ::   fsnd(:)            ! temporary arrays

REAL(r8), POINTER :: pfl_flx_total_col(:,:)
REAL(r8), POINTER :: pfl_flx_total_gcell(:,:)
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

   
CALL c2g(begc, endc, begl, endl, begg, endg, nlevsoi, pfl_flx_total_col, pfl_flx_total_gcell, &
         c2l_scale_type= 'unity', l2g_scale_type='unity')
   
      
ALLOCATE ( fsnd(begg:endg), stat=nerror)
   IF (nerror /= 0) THEN
      CALL prism_abort_proto( ncomp_id, 'send_fld_2pf', 'Failure in allocating fsnd' )
      RETURN
   END IF

! zero on unmasked points

   isec = dtime * ( get_nstep() -1 )

   DO jn = 1, vsnd
     fsnd = -999999._r8
     DO g1 = begg, endg
        fsnd(g1) = pfl_flx_total_gcell(g1,jn)*3.6_r8/dz(g1,jn)
     END DO
     !CMS the first 100th snd fields are allocated for COSMO-CLM coupling
     IF( ssnd(100+jn)%laction ) CALL oas_clm_snd( 100+jn, isec,fsnd(:),begg,endg, info ) 
   END DO
    
   isec = dtime * get_nstep()

   DEALLOCATE(fsnd)

CALL MPI_Barrier(kl_comm, nerror)

!CPS WRITE(6,*) 'oas_snd_pfl complete ...'

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE send_fld_2pfl

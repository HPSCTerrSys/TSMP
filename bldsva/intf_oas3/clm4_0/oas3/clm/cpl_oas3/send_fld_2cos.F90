SUBROUTINE send_fld_2cos(nstep, dtime)

!---------------------------------------------------------------------
! Description:
!  This routine sends the fluxes from CESM to COSMO
!
! Current Code Owner: TR32, Z4: Prabhakar Shrestha
!    phone: 0228733453
!    email: pshrestha@uni-bonn.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1.1        2011/11/28 Prabhakar Shrestha 
!   Modfied and Implemented in CLM3.5, Initial release
! 1.2.1        2012/09/18 Markus Uebel 
!   Inclusion of CO2 coupling (photosynthesis rate)
! 1.3.1        2013/02/01 Prabhakar Shrestha
!   nee used for CO2, albd and albi allocated for direct/diffuse albedos
!   Included aerodynamic resistance and surface temperature/moisture
! 1.3.1        2013/07/23 Fabian Gasper
!   Bug fix in albt_gcell for sending with multiple threads
!   This gives 2 options for coupling COSMO and CLM i.e either flux or transfer coefficients
! 2.1.0        2016/02/29 Prabhakar Shrestha
! Implementation for CESM 1.2.1
! @VERSION@    @DATE@     <Your name>
!  <Modification comments>         

!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! Declarations:
!
! Modules used:
 
USE oas_clm_vardef
!

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Local Variables:
INTEGER                          :: g              ! indices
INTEGER,  INTENT(IN)             :: nstep, dtime
INTEGER                          :: isec
! processor bounds indices
INTEGER, PARAMETER               ::   jps_co2fl  = 1 !CPScesm  8    !  net CO2 flux (now only photosynthesis rate) (umol CO2 m-2s-1)
INTEGER                          :: rank,info,jj, ji
REAL(KIND=r8), ALLOCATABLE       :: fsnd(:)      ! temporary arrays
!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine send_fld2cos 
!------------------------------------------------------------------------------

   CALL MPI_Comm_Rank(kl_comm, rank, nerror)  !CPScesm 
   IF (nerror /= 0) CALL prism_abort_proto(ncomp_id, 'MPI_Comm_Rank', 'Failure in send_fld_2cos') !CPScesm

   ALLOCATE ( fsnd(start1d:start1d+length1d-1), stat=nerror)  
   IF (nerror /= 0) THEN 
     CALL prism_abort_proto( ncomp_id, 'send_fld_2cos', 'Failure in allocating fsnd' )
     RETURN
   ENDIF

 ! zero on unmasked points
 fsnd = 1._r8*nstep
  
 isec = nstep*dtime
 PRINT*, "sendfld2cos: rank,isec, sending ...",rank, isec, MINVAL(fsnd), MAXVAL(fsnd)
 IF( ssnd(jps_co2fl)%laction )  CALL oas_clm_snd( jps_co2fl, isec, fsnd,start1d,start1d+length1d-1,info )

 CALL MPI_Barrier(kl_comm, nerror)

 DEALLOCATE(fsnd)
!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------


END SUBROUTINE send_fld_2cos

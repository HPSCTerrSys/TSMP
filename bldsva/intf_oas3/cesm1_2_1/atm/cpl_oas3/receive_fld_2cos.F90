SUBROUTINE receive_fld_2cos(nstep, dtime, lcoupled)

!---------------------------------------------------------------------
! Description:
!  This routine receives the atmospheric fields from COSMO and updates
!  the DATM a2x arrays 
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
! 1.2.1        2012/09/26 Markus Uebel, Prabhakar Shrestha 
!   CO2 coupling included
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
USE shr_kind_mod ,            ONLY : R8 => SHR_KIND_R8
!==============================================================================

IMPLICIT NONE

!==============================================================================

! Parameter List
INTEGER, INTENT(IN)                    ::   nstep                        ! time step
INTEGER, INTENT(IN)                    ::   dtime                        ! dt
LOGICAL, INTENT(OUT)                   ::  lcoupled                      ! CPS flag for coupling 
!
INTEGER, PARAMETER ::   jps_t   =  1            ! temperature
INTEGER, PARAMETER ::   jps_u   =  2            ! u wind
INTEGER, PARAMETER ::   jps_v   =  3            ! v wind
INTEGER, PARAMETER ::   jps_q   =  4            ! specific water vapor content
INTEGER, PARAMETER ::   jps_th  =  5            ! thickness of lowest level (m)
INTEGER, PARAMETER ::   jps_pr  =  6            ! surface pressure (Pa)
INTEGER, PARAMETER ::   jps_rs  =  7            ! direct shortwave downward radiation (W/m2)
INTEGER, PARAMETER ::   jps_fs  =  8            ! diffuse shortwave downward radiation (W/m2)
INTEGER, PARAMETER ::   jps_lw  =  9            ! longwave downward radiation (W/m2) 
INTEGER, PARAMETER ::   jps_cr  = 10            ! convective rain precipitation      (kg/m2*s)
INTEGER, PARAMETER ::   jps_cs  = 11            ! convective snow precipitation      (kg/m2*s)
INTEGER, PARAMETER ::   jps_gr  = 12            ! gridscale rain precipitation
INTEGER, PARAMETER ::   jps_gs  = 13            ! gridscale snow precipitation
INTEGER, PARAMETER ::   jps_gg  = 14            ! gridscale graupel precipitation
INTEGER, PARAMETER ::   jps_cp  = 15            ! total convective precipitation
INTEGER, PARAMETER ::   jps_gp  = 16            ! total gridscale precipitation
!MU (13.09.2012)
INTEGER, PARAMETER ::   jps_co2 = 17            ! CO2 partial pressure (Pa)
!MU (13.09.2012)

! Local Variables
INTEGER                                :: jn, isec, ier , begg,endg
INTEGER, DIMENSION(krcv)               :: nrcvinfo           ! OASIS info argument
INTEGER                                :: rank
REAL(KIND=r8), ALLOCATABLE             :: ztmp1(:)
!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine receive_fld2cos 
!------------------------------------------------------------------------------

 CALL MPI_Comm_Rank(kl_comm, rank, nerror)  !CPScesm 
 IF (nerror /= 0) CALL prism_abort_proto(ncomp_id, 'MPI_Comm_Rank', 'Failure in receive_fld_2cos') !CPScesm

 isec = dtime * nstep 

 ! Default value (changed when coupling)
 lcoupled = .FALSE.
 nrcvinfo = OASIS_idle
 begg     = start1d
 endg     = start1d+length1d-1

 ALLOCATE ( ztmp1(begg:endg), stat=nerror)
 IF (nerror /= 0) THEN
   CALL prism_abort_proto( ncomp_id, 'receive_fld_2cos', 'Failure in allocating ztmp1' )
   RETURN
 ENDIF
 
! Update DATM variable with the received input from COSMO
 DO jn = 1, krcv
   ztmp1 = -1._r8
   nrcvinfo = OASIS_idle
   IF(srcv(jn)%laction) CALL oas_clm_rcv( jn, isec, ztmp1,begg,endg, nrcvinfo(jn) )
   IF( nrcvinfo(jn)==OASIS_Rcv) THEN   
     !atm_a2l%forc_t(g) = ztmp1(g,1) + SHR_CONST_TKFRZ
     lcoupled = .TRUE.
   ENDIF
 ENDDO

 DEALLOCATE(ztmp1)

 !  Wait coupling result (involved PE and non involved PE)
 !CALL MPI_Bcast( lcoupled, 1, dummy, 0, kl_comm, ier )

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE receive_fld_2cos

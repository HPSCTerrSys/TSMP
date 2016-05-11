SUBROUTINE oas_clm_snd( kid, kstep, pdata,begg, endg, kinfo )

!---------------------------------------------------------------------
! Description:
!  This routine sends CESM fields to OASIS3 coupler at each coupling
!  time step defined in namcouple
!
! Current Code Owner: TR32, Z4: Prabhakar Shrestha
!    phone: 0228733453
!    email: pshrestha@uni-bonn.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
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

USE shr_kind_mod ,       ONLY : r8 => shr_kind_r8
USE oas_clm_vardef
USE mod_prism_put_proto

!==============================================================================

IMPLICIT NONE

!==============================================================================

! * Arguments
!
INTEGER,                 INTENT(IN)    :: kid    ! variable index in the array
INTEGER,                 INTENT(OUT)   :: kinfo  ! OASIS4 info argument
INTEGER,                 INTENT(IN)    :: kstep  ! CLM time-step in seconds
INTEGER,                 INTENT(IN)    :: begg,endg
REAL(KIND=r8), DIMENSION(begg:endg), INTENT(IN)     :: pdata
!
INTEGER                                             :: NULOUT

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine oas_clm_snd 
!------------------------------------------------------------------------------

 NULOUT=6
 !
 ! Call OASIS at each time step but field sent to other model only at coupling time step
 ! (accumulation otherwise, if asked in the SMIOC configuration file)
 !
 CALL prism_put_proto( ssnd(kid)%nid, kstep, pdata, kinfo )
 IF ( kinfo .NE. PRISM_Ok .AND. kinfo .LT. PRISM_Sent ) &
   CALL prism_abort_proto(kl_comm, 'oasclm', 'Failure in oas_clm_snd')
 IF ( IOASISDEBUGLVL == 1 )  &
  WRITE(NULOUT,*) "oas_clm_snd : ", kstep, ssnd(kid)%clname, MINVAL(pdata), MAXVAL(pdata)
!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE oas_clm_snd


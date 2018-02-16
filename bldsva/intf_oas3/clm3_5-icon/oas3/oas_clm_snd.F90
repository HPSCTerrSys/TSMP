SUBROUTINE oas_clm_snd( kid, kstep, pdata, kinfo )

!---------------------------------------------------------------------
! Description:
!  This routine sends CLM3.5 fields to OASIS3 coupler at each coupling
!  time step defined in namcouple
!
! References:
!  CEREFACS/ETH: E. Maisonnave, Edoward Davin
!
! Current Code Owner: TR32, Z4: Prabhakar Shrestha
!    phone: 0228733453
!    email: pshrestha@uni-bonn.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        2011/11/28 Prabhakar Shrestha 
!   Modfied and Implemented in CLM3.5, Initial release
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
REAL(KIND=r8), DIMENSION(ndlon,ndlat), INTENT(IN)   :: pdata
!
INTEGER                                             :: NULOUT
REAL(KIND=r8), DIMENSION(ndlon*ndlat)               :: ztmp1

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
 ztmp1(:) = RESHAPE (pdata(:,:), (/ndlon*ndlat/))
 CALL prism_put_proto( ssnd(kid)%nid, kstep, ztmp1, kinfo )

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE oas_clm_snd


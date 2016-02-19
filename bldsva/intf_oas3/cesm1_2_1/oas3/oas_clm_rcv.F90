SUBROUTINE oas_clm_rcv( kid, kstep, pdata, kinfo )

!---------------------------------------------------------------------
! Description:
!  This routine receiveds atmospheric field from OASIS3 coupler at
!  each timestep
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

USE shr_kind_mod , ONLY : r8 => shr_kind_r8
USE oas_clm_vardef
USE mod_prism_get_proto

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Arguments
INTEGER,                          INTENT(IN)        :: kid    ! variable intex in the array
INTEGER,                          INTENT(IN)        :: kstep  ! time-step in seconds
REAL(KIND=r8), DIMENSION(ndlon,ndlat), INTENT(OUT)  :: pdata
INTEGER,                          INTENT(OUT)       :: kinfo  ! OASIS info argument

! Local Variables
LOGICAL                                           :: llaction
INTEGER                                           :: NULOUT=6
REAL(r8), DIMENSION(ndlon*ndlat)                  :: ztmp1

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine oas_clm_rcv 
!------------------------------------------------------------------------------

 !  Masked point are not modified by OASIS
 !  = buffer set to 0 before calling oasis
 ztmp1=0._r8
 !
 ! receive data from OASIS
 !
 CALL prism_get_proto( srcv(kid)%nid, kstep, ztmp1, kinfo )         

 llaction = .false.
 IF( kinfo == PRISM_Recvd   .OR. kinfo == PRISM_FromRest .OR.   &
     kinfo == PRISM_RecvOut .OR. kinfo == PRISM_FromRestOut )   llaction = .TRUE.

 IF ( IOASISDEBUGLVL == 1 )  &
    WRITE(NULOUT,*) "oasclm: oas_clm_rcv: llaction, kinfo, kstep, clname: " , llaction, kinfo, kstep, srcv(kid)%clname

 IF ( llaction ) THEN

  ! Declare to calling routine that OASIS provided coupling field
  kinfo = OASIS_Rcv

 ! Update array which contains coupling field (only on valid shape)
   pdata(:,:) = RESHAPE (ztmp1(:), (/ndlon, ndlat/))
         
 ELSE
  ! Declare to calling routine that OASIS did not provide coupling field
  kinfo = OASIS_idle     
 ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------


END SUBROUTINE oas_clm_rcv

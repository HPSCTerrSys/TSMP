SUBROUTINE oas_cos_rcv( kid, kstep, pdata, kinfo )

!---------------------------------------------------------------------
! Description:
!  This routine call fields from OASIS3 coupler at each coupling time step.
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
!   Modfied and Implemented in COSMO4.11, Initial release
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

USE data_parameters,     ONLY:  ireals, iintegers
USE data_modelconfig, ONLY :  ie, je
USE oas_cos_vardef
USE mod_prism_get_proto

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Arguments

INTEGER(KIND=iintegers), INTENT(IN)   :: kid    ! variable intex in the array
INTEGER(KIND=iintegers), INTENT(IN)   :: kstep  ! ocean time-step in seconds
REAL(KIND=ireals), DIMENSION(ie,je), INTENT(OUT)   :: pdata
INTEGER(KIND=iintegers), INTENT(OUT)  :: kinfo  ! OASIS4 info argument
!
LOGICAL                   :: llaction
INTEGER(KIND=iintegers)   :: NULOUT=6

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine oas_cos_rcv 
!------------------------------------------------------------------------------

 !  Masked point are not modified by OASIS4
 !  = buffer set to 0 before calling oasis
 exfld=0._ireals
 !
 IF (lpe_cpl) THEN
   !
   ! receive local data from OASIS3 on every process involved in the coupling
   !
   CALL prism_get_proto( srcv(kid)%nid, kstep, exfld(nldi:nlei, nldj:nlej), kinfo )         

   llaction = .false.
   IF( kinfo == PRISM_Recvd   .OR. kinfo == PRISM_FromRest .OR.   &
      kinfo == PRISM_RecvOut .OR. kinfo == PRISM_FromRestOut )   llaction = .TRUE.

   IF ( IOASISDEBUGLVL == 1 )  WRITE(NULOUT,*) "oascos: oas_cos_rcv:  &
     llaction, kinfo, kstep, clname: " , llaction, kinfo, kstep, srcv(kid)%clname

   ! If coupling time step
   IF ( llaction ) THEN

     ! Declare to calling routine that OASIS provided coupling field
     kinfo = OASIS_Rcv

     ! Update array which contains coupling field (only on valid shape)
     pdata(nldi:nlei, nldj:nlej) = exfld(nldi:nlei, nldj:nlej)

   IF ( IOASISDEBUGLVL == 1 ) THEN
     WRITE(NULOUT,*) '****************'
     WRITE(NULOUT,*) 'oascos: prism_get: Incoming ', srcv(kid)%clname
     WRITE(NULOUT,*) 'oascos: prism_get: ivarid '  , srcv(kid)%nid
     WRITE(NULOUT,*) 'oascos: prism_get:   kstep', kstep
     WRITE(NULOUT,*) 'oascos: prism_get:   info ', kinfo
     WRITE(NULOUT,*) '     - Minimum value is ', MINVAL(pdata(nldi:nlei, nldj:nlej))
     WRITE(NULOUT,*) '     - Maximum value is ', MAXVAL(pdata(nldi:nlei, nldj:nlej))
     WRITE(NULOUT,*) '     -     Sum value is ', SUM(pdata(nldi:nlei, nldj:nlej))
     WRITE(NULOUT,*) '****************'
   ENDIF
         
   ELSE
     ! Declare to calling routine that OASIS did not provide coupling field
     kinfo = OASIS_idle     
   ENDIF

 ENDIF             ! lpe_cpl

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE oas_cos_rcv

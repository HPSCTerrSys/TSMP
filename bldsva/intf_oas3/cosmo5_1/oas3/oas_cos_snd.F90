SUBROUTINE oas_cos_snd( kid, kstep, pdata, kinfo )

!---------------------------------------------------------------------
! Description:
!  This routine sends fields to OASIS3 coupler at each coupling time step.
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

USE data_modelconfig, ONLY :  ie, je
USE oas_cos_vardef
USE mod_prism_put_proto

!==============================================================================

IMPLICIT NONE

!==============================================================================

! * Arguments
!
INTEGER(KIND=iintegers), INTENT(IN)    :: kid    ! variable intex in the array
INTEGER(KIND=iintegers), INTENT(OUT)   :: kinfo  ! OASIS info argument
INTEGER(KIND=iintegers), INTENT(IN)    :: kstep  ! ocean time-step in seconds
REAL(KIND=wp), DIMENSION(ie,je), INTENT(IN)   :: pdata
!
INTEGER(kind=iintegers) :: NULOUT=6

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine oas_cos_snd 
!------------------------------------------------------------------------------

 IF ( lpe_cpl ) THEN

   ! prepare array (only valid shape, without halos) for OASIS
   exfld(nldi:nlei, nldj:nlej) = pdata(nldi:nlei, nldj:nlej)

   ! Call OASIS at each time step but field sent to other model only at coupling time step
   ! (accumulation otherwise, if asked in the namcouple configuration file)
   !
   CALL prism_put_proto( ssnd(kid)%nid, kstep, exfld(nldi:nlei, nldj:nlej), kinfo )

   IF ( IOASISDEBUGLVL == 1 ) THEN
     WRITE(NULOUT,*) '****************'
     WRITE(NULOUT,*) 'oascos: prism_put: Incoming ', ssnd(kid)%clname
     WRITE(NULOUT,*) 'oascos: prism_put: ivarid '  , ssnd(kid)%nid
     WRITE(NULOUT,*) 'oascos: prism_put:   kstep', kstep
     WRITE(NULOUT,*) 'oascos: prism_put:   info ', kinfo
     WRITE(NULOUT,*) '     - Minimum value is ', MINVAL(pdata(nldi:nlei, nldj:nlej))
     WRITE(NULOUT,*) '     - Maximum value is ', MAXVAL(pdata(nldi:nlei, nldj:nlej))
     WRITE(NULOUT,*) '     -     Sum value is ', SUM(pdata(nldi:nlei, nldj:nlej))
     WRITE(NULOUT,*) '****************'
   ENDIF

 ENDIF       ! lpe_cpl

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------


END SUBROUTINE oas_cos_snd


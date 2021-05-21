! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
! Atmospheric models
!  1 : L-band, Pellarin
!  2 : Ulaby
!===========================================================================

SUBROUTINE ATMPELLARIN

! Purpose :
! -------
!   Calculate atmospheric opacity and up- and downwelling 
!   atmospheric radiation at L-band.
!   For frequencies above 10 GHz, surface water vapor density should 
!   be accounted for
   

! Modifications :
! -------------
!   24-Sept-2006 Thomas Holmes, ECMWF
!   January 2008 Patricia de Rosnay, ECMWF
! End Modifications

! Internal variables :
!  GOSSAT : ATMOSPHERIC LOSS FACTOR: GOSSAT = EXP( -tau_atm / costheta)
!  TAEQ : equivalent layer temperature [K]
!---------------------------------------------------------------------------

USE PARKIND1, ONLY : JPRM
USE YOMCMEMPAR, ONLY : costheta, LGPRINT,LOFIELDEXP
USE YOMCMEMATM, ONLY : tair, Z,  t_sky, tau_atm, tb_ad, tb_au

IMPLICIT NONE

REAL(KIND=JPRM) :: GOSSAT, TAEQ
!---------------------------------------------------------------------------

! 1. Zenith atmospheric opacity

tau_atm = EXP( -3.926 - 0.2211 * Z - 0.00369 *tair)


! 2. Calculate up- and downward atmospheric radiation

GOSSAT = EXP(-tau_atm/ costheta)


TAEQ = EXP( 4.927 + 0.002195 * tair)

tb_ad = TAEQ*(1. - GOSSAT) + t_sky * GOSSAT


tb_au =  TAEQ*(1. - GOSSAT) 

SELECT CASE (LOFIELDEXP)

     CASE (.TRUE.)
      tb_au = 0.0
ENDSELECT

END SUBROUTINE ATMPELLARIN


!===========================================================================

SUBROUTINE ATMULABY

! Purpose :
! -------
!   Calculate atmospheric opacity and up- and downwelling atmospheric 
!   radiation for L-band
   
!   24-Sept-2006 Thomas Holmes   *ECMWF*

! LOSSAT=ATMOSPHERIC LOSS FACTOR: LOSSAT(DB)=TAUAT(DB)*SEC(TETA0)
!---------------------------------------------------------------------------

USE PARKIND1, ONLY : JPRM
USE YOMCMEMPAR, ONLY : costheta, LGPRINT
USE YOMCMEMATM, ONLY : tair, t_sky, tau_atm, tb_ad, tb_au

IMPLICIT NONE

REAL(KIND=JPRM) :: LOSSAT
!---------------------------------------------------------------------------

! 1. Zenith atmospheric opacity

CALL ULABY_TAU_ATM

! 2. Calculate up- and downward atmospheric radiation

LOSSAT = EXP((tau_atm/costheta)/4.34)

tb_ad = tair*(1. -1./LOSSAT) + t_sky/LOSSAT

tb_au = tair*(1. -1./LOSSAT) 

END SUBROUTINE ATMULABY

!===========================================================================

SUBROUTINE ULABY_TAU_ATM

! Purpose :
! -------
!   Calculate zenith atmospheric opacity
   
! Modifications :
! -------------
!   24-Sept-2006 Thomas Holmes   *ECMWF*
! End Modifications

!  ROSWAT=SURFACE-WATER-VAPOR-DENSITY (ROSWAT=7.5G.M-3)
!---------------------------------------------------------------------------

USE PARKIND1, ONLY : JPRM
USE YOMCMEMPAR, ONLY : fghz, LGPRINT
USE YOMCMEMATM, ONLY : tau_atm

IMPLICIT NONE

REAL(KIND=JPRM), PARAMETER :: ROSWAT = 7.5
REAL(KIND=JPRM) :: X
!---------------------------------------------------------------------------

tau_atm = 0.

X=ABS((fghz-90.)/90.)
IF (X <= 0.1) THEN
    tau_atm=0.17+0.06*ROSWAT
ENDIF
X=ABS((fghz-35.)/35.)
IF (X <= 0.1) THEN
    tau_atm=0.17+0.013*ROSWAT
ENDIF
X=ABS((fghz-22.235)/22.235)
IF (X <= 0.1) THEN
    tau_atm=0.11+0.048*ROSWAT
ENDIF
X=ABS((fghz-15.)/15.)
IF (X <= 0.1) THEN
    tau_atm=0.055+0.004*ROSWAT
ENDIF
IF (fghz <= 11.) THEN
    tau_atm=0.04
ENDIF

IF (tau_atm == 0.) THEN
    CALL ABOR1('tau_atm == 0')
ENDIF

END SUBROUTINE ULABY_TAU_ATM


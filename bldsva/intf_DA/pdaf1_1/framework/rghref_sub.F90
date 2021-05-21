! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
!---------------------------------------------------------------------------
! Different models to calculate Rough Surface Emissivity
!  1,2,4,5) 'Choudhury', 'Wsimple', 'Wtexture', 'Wigneron' use RGHCHOU 
!  3) 'Wegmueller' uses RGHWEGM 
!---------------------------------------------------------------------------

SUBROUTINE RGHCHOU

! Purpose :
!   Calculate Rough Surface Emissivity

  
!---------------------------------------------------------------------------

USE YOMCMEMPAR, ONLY : costheta, LGPRINT, ip_Q
USE YOMCMEMSOIL, ONLY : r_r, r_s, Nrv,Nrh,hrmodel 
USE YOMCMEMFIELDS, ONLY : JJ, N

IMPLICIT NONE
!---------------------------------------------------------------------------

r_r(1) = (ip_Q * r_s(2) + (1.-ip_Q) * r_s(1)) * exp(-hrmodel * costheta**Nrh)
r_r(2) = (ip_Q * r_s(1) + (1.-ip_Q) * r_s(2)) * exp(-hrmodel * costheta**Nrv)



END SUBROUTINE RGHCHOU

!===========================================================================

SUBROUTINE RGHWEGM

! Purpose :
! calculate Rough Surface Emissivity based on smooth surface reflectivity (H only)
!---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRM

USE YOMCMEMPAR ,ONLY : theta, costheta, LGPRINT
USE YOMCMEMSOIL, ONLY : r_r , r_s, hrmodel 
USE YOMCMEMFIELDS, ONLY : JJ, N, fh

IMPLICIT NONE
!---------------------------------------------------------------------------

r_r(1) = r_s(1) * exp(-1.0 * hrmodel ** (sqrt(0.10 * costheta)))

IF (hrmodel == 0.) THEN
  r_r(2) = r_s(2) 
ELSE
  IF (theta <= 60.0) THEN
    r_r(2) = r_r(1) * (costheta ** 0.655)
  ELSEIF (theta <= 70.0) THEN
    r_r(2) = r_r(1) * (0.635_JPRM - 0.0014_JPRM*(theta-60.0_JPRM))
  ELSE
   CALL ABOR1('Incidence angle not in [0;70] degrees. not suitable for rough soil reflectivity')
  ENDIF
ENDIF

END SUBROUTINE RGHWEGM

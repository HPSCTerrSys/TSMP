! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
!
SUBROUTINE CMEM_ATM

! Purpose :
! -------
! Compute Atmospheric contribution to the microwave emission
   

! Internal variables
! ------------------
! Authors :
!     Patricia de Rosnay  Feb 2008, January 2009

!------------------------------------------------------------------------------

USE PARKIND1, ONLY : JPIM, JPRM
USE YOMLUN , ONLY : NULOUT

USE YOMCMEMPAR, ONLY : CIATM, LGPRINT
USE YOMCMEMFIELDS, ONLY: N, JJ, ftair, ftau_atm, ftb_au,ftb_ad
USE YOMCMEMATM, ONLY : fZ,tb_au,tb_ad,tau_atm,tair,Z

IMPLICIT NONE
INTEGER(KIND=JPIM) :: i


!------------------------------------------------------------------------------

  DO JJ = 1, N
    
    tau_atm = 0._JPRM
    tb_au = 0._JPRM
    tb_ad = 0._JPRM
    tair = ftair (JJ)

    SELECT CASE (CIATM)

      CASE ( 'Pellarin' )

        Z = fZ(JJ)
        CALL ATMPELLARIN

      CASE ( 'Ulaby' )

        CALL ATMULABY

    END SELECT

    ftau_atm(JJ) = tau_atm
    ftb_au(JJ) = tb_au
    ftb_ad(JJ) = tb_ad


  ENDDO


IF (LGPRINT) WRITE(NULOUT,*) 'End of CMEM Atmospheric Module'

       
END SUBROUTINE CMEM_ATM


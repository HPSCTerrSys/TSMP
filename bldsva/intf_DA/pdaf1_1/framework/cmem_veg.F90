! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
SUBROUTINE CMEM_VEG

! Purpose :
! -------
!  Vegetation model

! Author :
! ------
!   23-Aug-2006 Thomas Holmes   *ECMWF*
!   13 December 2013 Patricia de Rosnay (ECMWF) and Nathalie Gautier (Env. Canada)
!               account for frozen vegetation

! External :
!   dielsal_sub.F90
!   veg_sub.F90
!---------------------------------------------------------------------------
USE PARKIND1, ONLY : JPIM, JPRM
USE YOMLUN   , ONLY : NULOUT
USE YOMCMEMPAR, ONLY : CIVEG, tfreeze,LGPRINT
USE YOMCMEMVEG, ONLY : t_veg, tau_veg, tb_veg, sal_vw, eps_vw, w_eff

IMPLICIT NONE

INTEGER(KIND=JPIM) :: i
INTEGER(KIND=JPIM) :: medium, isal
!---------------------------------------------------------------------------

VEG: SELECT CASE (CIVEG)

  CASE ( 'No' )
    tau_veg(:) = (/0.,0./)
    tb_veg(:) = (/0.,0./)
    
  CASE ('Kirdyashev','Wegmueller','Wigneron','Jackson' ) 
    
    ! 2. Calculate vegetation opacity (tau_veg)
    !    -------------------------------------
    medium = 0
    isal = 0
    SELECT CASE (CIVEG)
    CASE ( 'Kirdyashev' )
      Frozen_veg_Ki: SELECT CASE ( (t_veg-tfreeze < -0.5) )
      CASE ( .TRUE. ) Frozen_veg_Ki
         CALL DIEL_ICE (t_veg,eps_vw)
      CASE DEFAULT  
         CALL DIEL_WAT ( medium, isal,(t_veg-tfreeze), sal_vw, eps_vw)
      END SELECT Frozen_veg_Ki
      CALL vegkird  
    CASE ( 'Wegmueller' )
      Frozen_veg_We: SELECT CASE ( (t_veg-tfreeze < -0.5) )
      CASE ( .TRUE. ) Frozen_veg_We
         CALL DIEL_ICE (t_veg,eps_vw)
      CASE DEFAULT
         CALL DIEL_WAT ( medium, isal, (t_veg-tfreeze), sal_vw, eps_vw)
      END SELECT Frozen_veg_We
      CALL vegwegm
    CASE ( 'Wigneron' )
      CALL vegwign
    CASE ( 'Jackson' )
      CALL vegjack
    ENDSELECT
    IF (LGPRINT) WRITE(NULOUT,*) 'Tau veg is:' ,tau_veg

    ! 3. Calculate microwave emission from vegetation (tb_veg)
    !    -----------------------------------------------------
    DO i = 1, 2
      tb_veg(i) = (1. - w_eff(i)) * (1. - exp(-tau_veg(i))) * t_veg
    ENDDO
    IF (LGPRINT) WRITE(NULOUT,*) '--- T veg :' ,t_veg

ENDSELECT VEG

END SUBROUTINE CMEM_VEG
